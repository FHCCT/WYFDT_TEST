#include "dt.h"
#include "dt_bw.h"

#include <stdexcept>

namespace dt {
int64_t& BWScratchTable::operator[](int key) {
    size_t cap = large.empty() ? small.size() : large.size();
    Slot* slots = large.empty() ? small.data() : large.data();
    auto position = [](int k, size_t mask) { return (uint32_t(k) * 2654435761u) & mask; };
    size_t pos = position(key, cap - 1);
    while (slots[pos].key != -1 && slots[pos].key != key) pos = (pos + 1) & (cap - 1);
    if (slots[pos].key == key) return slots[pos].value;
    if ((count + 1) * 2 > cap) {
        grow();
        return (*this)[key];
    }
    slots[pos].key = key; ++count;
    return slots[pos].value;
}

void BWScratchTable::grow() {
    const size_t cap = large.empty() ? small.size() : large.size();
    const Slot* slots = large.empty() ? small.data() : large.data();
    std::vector<Slot> grown(cap * 2);
    for (size_t i = 0; i < cap; ++i) if (slots[i].key != -1) {
        size_t dst = (uint32_t(slots[i].key) * 2654435761u) & (grown.size() - 1);
        while (grown[dst].key != -1) dst = (dst + 1) & (grown.size() - 1);
        grown[dst] = slots[i];
    }
    large.swap(grown);
}

void BWScratchTable::clear() {
    // Retain common capacities; only shrink exceptional cavities. This avoids
    // allocating again when successive ordinary cavities straddle a size limit.
    if (!large.empty() && large.size() > std::max(size_t(512), count * 4)) {
        size_t target = 256;
        while (target < count * 2) target *= 2;
        large.resize(target);
    }
    if (large.empty()) small.fill(Slot());
    else std::fill(large.begin(), large.end(), Slot());
    count = 0;
}

void BWEdgeTable::clear() {
    if (!large.empty() && large.size() > std::max(small.size(), count * 4)) {
        if (count <= small.size() / 2) large.clear();
        else {
            size_t target = small.size();
            while (target < count * 2) target *= 2;
            large.resize(target);
        }
    }
    if (large.empty()) small.fill(Slot());
    else std::fill(large.begin(), large.end(), Slot());
    count = 0;
}

BWEdgeTable::Slot& BWEdgeTable::insert(uint64_t key, int face, int side, bool& inserted) {
    auto position = [](uint64_t value, size_t mask) {
        value ^= value >> 33;
        value *= 0xff51afd7ed558ccdULL;
        value ^= value >> 33;
        return size_t(value) & mask;
    };
    size_t cap = large.empty() ? small.size() : large.size();
    Slot* slots = large.empty() ? small.data() : large.data();
    size_t pos = position(key, cap - 1);
    while (slots[pos].key != UINT64_MAX && slots[pos].key != key) pos = (pos + 1) & (cap - 1);
    if (slots[pos].key == key) { inserted = false; return slots[pos]; }
    if ((count + 1) * 2 > cap) {
        std::vector<Slot> grown(cap * 2);
        for (size_t i = 0; i < cap; ++i) if (slots[i].key != UINT64_MAX) {
            size_t dst = position(slots[i].key, grown.size() - 1);
            while (grown[dst].key != UINT64_MAX) dst = (dst + 1) & (grown.size() - 1);
            grown[dst] = slots[i];
        }
        large.swap(grown); cap = large.size(); slots = large.data();
        pos = position(key, cap - 1);
        while (slots[pos].key != UINT64_MAX) pos = (pos + 1) & (cap - 1);
    }
    Slot& slot = slots[pos];
    slot.key = key; slot.face = face; slot.side = static_cast<uint8_t>(side); slot.paired = false;
    ++count; inserted = true;
    return slot;
}

namespace {
constexpr int64_t cavityFlag = int64_t(1) << 29;
constexpr int64_t testedFlag = int64_t(1) << 30;
constexpr int64_t readFlag = 2;
constexpr int64_t writeFlag = 4;

bool validTet(DT& m, int t) {
    return t >= 0 && t < static_cast<int>(m.Elems.size()) && m.Elems[t].info >= 0;
}

int64_t* checkedTetInfo(DT& m, BWPlan& p, int t) {
    if (!validTet(m, t)) { p.status = BWStatus::Stale; return nullptr; }
    int64_t& flags = p.tetInfo[t];
    if (p.request.trackAccess && !(flags & readFlag)) { flags |= readFlag; p.readTets.push_back(t); }
    return &flags;
}

bool recordRead(DT& m, BWPlan& p, int t) {
    if (!validTet(m, t)) { p.status = BWStatus::Stale; return false; }
    if (!p.request.trackAccess) return true;
    int64_t& flags = p.tetInfo[t];
    if (!(flags & readFlag)) { flags |= readFlag; p.readTets.push_back(t); }
    return true;
}

void appendCavity(BWPlan& p, int t, int64_t& flags) {
    if (flags & cavityFlag) return;
    flags |= cavityFlag | (int64_t(p.workingCavity.size()) << 32);
    p.workingCavity.push_back(t);
    p.boundaryFaces.push_back(0);
}

void appendCavity(BWPlan& p, int t) {
    appendCavity(p, t, p.tetInfo[t]);
}

void appendWriteTet(BWPlan& p, int tet, int64_t& flags) {
    if (!(flags & writeFlag)) { flags |= writeFlag; p.writeTets.push_back(tet); }
}

void appendWriteNode(BWPlan& p, int node) {
    int64_t& flags = p.nodeInfo[node];
    if (!(flags & writeFlag)) { flags |= writeFlag; p.writeNodes.push_back(node); }
}

// An absent boundary vertex proves that this face is not in the frozen table.
// Positive membership is only a filter: it must still use the original hash.
bool isBoundaryTriangle(DT& m, const BWRequest& request, int a, int b, int c) {
    if (request.info == 2 && request.boundaryNodes) {
        const auto& nodes = *request.boundaryNodes;
        if (static_cast<size_t>(a) >= nodes.size() || !nodes[a] ||
            static_cast<size_t>(b) >= nodes.size() || !nodes[b] ||
            static_cast<size_t>(c) >= nodes.size() || !nodes[c]) return false;
    }
    return m.isBndTri(a, b, c);
}

// Location is read-only and runs against the frozen mesh. Walk-only cells do
// not affect the final plan once the containing cell has been determined.
int locateRequest(DT& m, BWPlan& p, int& tet) {
    if (!validTet(m, tet)) { p.status = BWStatus::Stale; return 100; }
    const int init = tet;
    if (m.ishulltet(tet)) tet = m.getNeig(tet, 3);
    int prev = -1, zero[2] = {-1, -1};
    for (int attempt = 0; attempt < 10000; ++attempt) {
        if (!validTet(m, tet)) { p.status = BWStatus::Stale; return 100; }
        if (m.ishulltet(tet) || (prev != -1 && tet == init)) return 1;
        double smallest = DBL_MAX;
        int nextFace = -1;
        zero[0] = zero[1] = -1;
        for (int f = 0; f < 4; ++f) {
            if (m.getNeig(tet, f) == prev) continue;
            int ia, ib, ic, id;
            DNC(f, ia, ib, ic, id);
            const int b = m.Elems[tet].form[ib], c = m.Elems[tet].form[ic], d = m.Elems[tet].form[id];
            const double ori = GEOM_FUNC::orient3d(m.Nodes[b].pt, m.Nodes[c].pt, m.Nodes[d].pt, p.request.point.data());
            if (ori < 0) {
                const double value = ori / std::max(m.calArea(m.Nodes[b].pt, m.Nodes[c].pt, m.Nodes[d].pt), 1e-30);
                if (value < smallest) { smallest = value; nextFace = f; }
            } else if (ori == 0) {
                if (zero[0] == -1) zero[0] = f;
                else if (zero[1] == -1) zero[1] = f;
                else {
                    if (ib != zero[0] && ib != zero[1]) p.duplicateNode = b;
                    else if (ic != zero[0] && ic != zero[1]) p.duplicateNode = c;
                    else p.duplicateNode = d;
                    p.status = BWStatus::Duplicate;
                    return -1;
                }
            }
        }
        if (nextFace != -1) { prev = tet; tet = m.getNeig(tet, nextFace); continue; }
        if (zero[0] != -1 && zero[1] == -1) return (zero[0] + 1) << 4;
        if (zero[0] != -1 && zero[1] != -1) {
            int b = -1, c = -1;
            for (int i = 3; i >= 0; --i) if (i != zero[0] && i != zero[1]) {
                if (b == -1) b = i; else c = i;
            }
            return (b << 2) | c;
        }
        return 1;
    }
    p.status = BWStatus::Rejected;
    return 100;
}

// The shell walk owns no locks and writes no mesh markers.
bool readShell(DT& m, BWPlan& p, int tet, int p2, int p3, std::vector<int>& shell) {
    shell.clear();
    if (p2 < 0 || p3 < 0 || p2 == p3 || !recordRead(m, p, tet)) return false;
    int p0 = 0, p1 = 0;
    DDNC(p0, p1, p2, p3);
    const int first = m.Elems[tet].form[p0];
    int end = m.Elems[tet].form[p1];
    const int a = m.Elems[tet].form[p2], b = m.Elems[tet].form[p3];
    shell.push_back(tet);
    while (end != first) {
        const int face = m.getNeigOrd(tet, p0);
        tet = m.getNeig(tet, p0);
        if (!recordRead(m, p, tet) || shell.size() > 10000) return false;
        end = m.Elems[tet].form[face];
        for (int i = 0; i < 4; ++i)
            if (i != face && m.Elems[tet].form[i] != a && m.Elems[tet].form[i] != b) { p0 = i; break; }
        shell.push_back(tet);
    }
    return true;
}

bool matchSize(DT& m, BWPlan& p, int node, int64_t& flags) {
    const BWRequest& r = p.request;
    if ((flags & 1) || node == m.ghost) return true;
    flags |= 1;
    const double distance = m.distance(p.request.point.data(), m.Nodes[node].pt);
    if (r.anisotropic) {
        if (m.minEdge != -1 && distance < m.minEdge * m.ani_lower) return false;
        double a[6], b[6];
        for (int i = 0; i < 6; ++i) { a[i] = r.metric[i]; b[i] = m.AniSol[node][i]; }
        return m.cal_ani_length(p.request.point.data(), m.Nodes[node].pt, a, b) >= m.ani_lower;
    }
    if (m.minEdge != -1 && distance < m.minEdge) return false;
    return !(distance < m.meshSize * r.space && distance < m.meshSize * m.Nodes[node].space);
}
}

BWRequest makeBWRequest(DT& m, int node, const std::vector<int>& seeds, int info) {
    BWRequest r;
    r.node = node; r.seeds = seeds; r.info = info; r.initShell = m.initBWshell;
    for (int i = 0; i < 3; ++i) r.point[i] = m.Nodes[node].pt[i];
    r.space = m.Nodes[node].space; r.boundary = m.isbndpnt(node);
    r.anisotropic = !m.AniSol.empty();
    if (r.anisotropic) r.metric = m.AniSol[node];
    return r;
}

BWStatus findBWCavity(DT& m, BWPlan& p) {
    p.status = BWStatus::Ready; p.duplicateNode = -1; p.ghostTet = -1;
    p.readTets.clear(); p.writeTets.clear(); p.writeNodes.clear(); p.cavity.clear();
    p.faces.clear(); p.newElements.clear(); p.workingCavity.clear(); p.boundaryFaces.clear(); p.lockNodes.clear();
    p.tetInfo.clear(); p.nodeInfo.clear();
    const BWRequest& r = p.request;
    if (r.seeds.empty() || (r.info == 0 && r.node < 0)) return p.status = BWStatus::Rejected;
    if (r.seeds.size() == 1) {
        int tet = r.seeds[0];
        const int loc = locateRequest(m, p, tet);
        if (p.status != BWStatus::Ready) return p.status;
        if (loc == 100) return p.status = BWStatus::Rejected;
        int64_t* initialInfo = checkedTetInfo(m, p, tet);
        if (!initialInfo) return p.status;
        if (loc == 1) appendCavity(p, tet, *initialInfo);
        else if (loc > 15) {
            appendCavity(p, tet, *initialInfo);
            const int next = m.getNeig(tet, (loc >> 4) - 1);
            int64_t* nextInfo = checkedTetInfo(m, p, next);
            if (!nextInfo) return p.status;
            appendCavity(p, next, *nextInfo);
        } else {
            auto& shell = p.shell;
            if (!readShell(m, p, tet, (loc & 12) >> 2, loc & 3, shell)) return p.status = BWStatus::Rejected;
            for (int t : shell) appendCavity(p, t);
        }
    } else for (int t : r.seeds) {
        int64_t* initialInfo = checkedTetInfo(m, p, t);
        if (!initialInfo) return p.status;
        appendCavity(p, t, *initialInfo);
    }

    // Request-wide policy is invariant during the entire cavity search.
    const bool checkSpacing = !r.boundary && (r.info == 2 || (r.info == 3 && r.anisotropic));
    for (size_t head = 0; head < p.workingCavity.size(); ++head) {
        const int tet = p.workingCavity[head];
        // Every queue entry was already validated and recorded when enqueued.
        // The mesh stays frozen until all plans in this batch are complete.
        for (int f = 0; f < 4; ++f) {
            const int node = m.Elems[tet].form[f];
            int64_t& nodeFlags = p.nodeInfo[node];
            if (!(nodeFlags & 2)) { nodeFlags |= 2; p.lockNodes.push_back(node); }
            if (checkSpacing && !matchSize(m, p, node, nodeFlags)) return p.status = BWStatus::Rejected;
            const int next = m.getNeig(tet, f);
            int64_t* nextInfo = checkedTetInfo(m, p, next);
            if (!nextInfo) return p.status;
            if (r.info != 0) {
                int ia, ib, ic, id;
                DFC(f, ia, ib, ic, id);
                const int b = m.Elems[tet].form[ib], c = m.Elems[tet].form[ic], d = m.Elems[tet].form[id];
                bool constrained = isBoundaryTriangle(m, r, b, c, d);
                if (!constrained && r.info == 3) {
                    auto singleEdge = [&](int u, int v) {
                        const int* edge = m.BndEdg.find(u, v);
                        return edge && m.SurEdgs[*edge].constrain > 0 && m.SurEdgs[*edge].face.size() == 1;
                    };
                    constrained = singleEdge(b, c) || singleEdge(b, d) || singleEdge(c, d);
                }
                if (constrained) { *nextInfo |= testedFlag; p.boundaryFaces[head] |= 1 << f; continue; }
            }
            if (*nextInfo & cavityFlag) continue;
            bool include = false;
            bool infoMayHaveMoved = false;
            if (!(*nextInfo & testedFlag) && !(r.info != 0 && m.ishulltet(next))) {
                // A hull second-neighbor lookup can grow tetInfo. Reacquire the
                // reference afterward; ordinary finite predicates need no lookup.
                *nextInfo |= testedFlag;
                const int a = m.Elems[next].form[0], b = m.Elems[next].form[1];
                const int c = m.Elems[next].form[2], d = m.Elems[next].form[3];
                if (!m.ishulltet(next)) {
                    const double value = r.info == 0 ? m.insphere_s(a, b, c, d, r.node)
                        : GEOM_FUNC::insphere(m.Nodes[a].pt, m.Nodes[b].pt, m.Nodes[c].pt, m.Nodes[d].pt, p.request.point.data());
                    include = value <= 0;
                } else {
                    const double value = GEOM_FUNC::orient3d(m.Nodes[a].pt, m.Nodes[b].pt, m.Nodes[c].pt, p.request.point.data());
                    include = value < 0;
                    if (value == 0) {
                        const int inner = m.getNeig(next, 3);
                        if (!recordRead(m, p, inner)) return p.status;
                        infoMayHaveMoved = true;
                        const int* v = m.Elems[inner].form;
                        include = m.insphere_s(v[0], v[1], v[2], v[3], r.node) < 0;
                    }
                }
            }
            if (include) {
                if (infoMayHaveMoved) nextInfo = &p.tetInfo[next];
                appendCavity(p, next, *nextInfo);
            }
            else p.boundaryFaces[head] |= 1 << f;
        }
    }
    return p.status;
}

BWStatus adjustBWCavity(DT& m, BWPlan& p) {
    if (p.status != BWStatus::Ready) return p.status;
    int ia, ib, ic, id, oldtet, outtet, outord, k;
	auto& inQ = p.adjustQueued;
	auto& q = p.adjustQueue;
	/************************** fix the cavity C(p) ***********************/
	while (1) {
		inQ.assign(p.workingCavity.size(), 0);
		q.clear();

		for (int i = 0; i < p.workingCavity.size(); ++i) {
			if (p.workingCavity[i] != -1) {
				q.push_back(i);
				inQ[i] = 1;
			}
		}

		while (!q.empty()) {
			int i = q.back(); q.pop_back();
			inQ[i] = 0;

			if (p.workingCavity[i] == -1) continue;      // already deleted

			int oldtet = p.workingCavity[i];

			for (int j = 0; j < 4; ++j) {
				if (!m.get_bit(p.boundaryFaces[i], j)) continue;

				int outtet = m.getNeig(oldtet, j);
                if (!recordRead(m, p, outtet)) return p.status;

				DNC(j, ia, ib, ic, id);
				int b = m.Elems[oldtet].form[ib];
				int c = m.Elems[oldtet].form[ic];
				int d = m.Elems[oldtet].form[id];

				if (m.get_bit(p.tetInfo[outtet], 29)) {
					outord = m.getNeigOrd(oldtet, j);
					if (!isBoundaryTriangle(m, p.request, b, c, d)) {
						m.clear_bit(p.boundaryFaces[p.tetInfo[outtet] >> 32], outord);
						m.clear_bit(p.boundaryFaces[i], j);
						continue;
					}
				}

				if (b == m.ghost || c == m.ghost || d == m.ghost)
					continue;

				double ori = dt::GEOM_FUNC::orient3d(p.request.point.data(), m.Nodes[b].pt, m.Nodes[d].pt, m.Nodes[c].pt);
				if (ori <= m.minVolume_bw) {
					for (int k = 0; k < 4; ++k) {
						outtet = m.getNeig(oldtet, k);
                            if (outtet != -1 && !recordRead(m, p, outtet)) return p.status;
						if (outtet == -1) continue;
						if (!m.get_bit(p.tetInfo[outtet], 29)) {
							p.tetInfo[outtet] &= 0x1FFFFFFF;
						}
						else {
							outord = m.getNeigOrd(oldtet, k);
							int ni = (int)(p.tetInfo[outtet] >> 32);
							m.set_bit(p.boundaryFaces[ni], outord);

							if (!inQ[ni]) { q.push_back(ni); inQ[ni] = 1; }
						}
					}

					p.workingCavity[(int)(p.tetInfo[oldtet] >> 32)] = -1;
					p.tetInfo[oldtet] &= 0x1FFFFFFF;

					break;
				}
			}
		}

		bool skip = true;
		if (p.request.info == 3) {
			static constexpr uint8_t faceKill[4] = {
				(1u << 3) | (1u << 4) | (1u << 5), // face 0
				(1u << 1) | (1u << 2) | (1u << 5), // face 1
				(1u << 0) | (1u << 2) | (1u << 4), // face 2
				(1u << 0) | (1u << 1) | (1u << 3)  // face 3
			};
			auto EdgeKey = [](uint32_t a, uint32_t b) -> uint64_t {
				if (a > b) std::swap(a, b);
				return (uint64_t)a << 32 | (uint64_t)b;
				};

			std::unordered_map<uint64_t, int> candEdges;
			for (int i = 0; i < p.workingCavity.size(); ++i) {
				if (p.workingCavity[i] == -1) continue;
				int oldtet = p.workingCavity[i];

				uint8_t edgeMask = 0x3F; // 6 edges all possible
				for (int f = 0; f < 4; ++f) {
					if (m.get_bit(p.boundaryFaces[i], f)) edgeMask &= (uint8_t)~faceKill[f];
				}

				for (int e = 0; e < 6; ++e) {
					int u = m.Elems[oldtet].form[Egid[e][0]];
					int v = m.Elems[oldtet].form[Egid[e][1]];
					uint64_t key = EdgeKey(u, v);

					auto it = candEdges.find(key);
					if (it != candEdges.end() && it->second == -1)
						continue;

					if (((edgeMask >> e) & 1u) == 0u) {// Case 1: excluded by boundary faces -> kill forever
						candEdges[key] = -1;
						continue;
					}

					if (!m.isBndEdg(u, v)) {  // Case 2: not a recovered/boundary edge -> kill forever
						candEdges[key] = -1;
						continue;
					}

					if (it == candEdges.end()) {
						candEdges.emplace(key, oldtet);
					}
				}
			}

			for (auto it : candEdges) {
				if (it.second == -1) continue;
				uint64_t key = it.first;
				uint32_t u = static_cast<uint32_t>(key >> 32);
				uint32_t v = static_cast<uint32_t>(key & 0xFFFFFFFFu);

				oldtet = it.second;
				auto& shell = p.shell;

				if (!readShell(m, p, oldtet, m.isNod_in_Tet(u, oldtet), m.isNod_in_Tet(v, oldtet), shell)) return p.status = BWStatus::Rejected;
				int numIn = 0;
				for (int k = 0; k < shell.size(); k++) {
					if (m.get_bit(p.tetInfo[shell[k]], 29)) {
						numIn++;
					}
				}
				if (numIn == shell.size()) {
					for (int t = p.workingCavity.size() - 1; t >= 0; t--) {
						oldtet = p.workingCavity[t];
						if (std::find(shell.begin(), shell.end(), oldtet) != shell.end()) {
							for (k = 0; k < 4; k++) {
								outtet = m.getNeig(oldtet, k);
                            if (outtet != -1 && !recordRead(m, p, outtet)) return p.status;
								if (outtet == -1)
									continue;
								if (!m.get_bit(p.tetInfo[outtet], 29)) {
									p.tetInfo[outtet] &= 0x1FFFFFFF;
								}
								else {
									outord = m.getNeigOrd(oldtet, k);
									m.set_bit(p.boundaryFaces[p.tetInfo[outtet] >> 32], outord);
								}
							}
							p.workingCavity[t] = -1;
							p.tetInfo[oldtet] &= 0x1FFFFFFF;
							skip = false;//try to fix again
							break;
						}
					}
				}
			}
		}

		if (skip)
			break;
	}

    return p.status;
}


bool pairBWBoundary(BWPlan& p) {
    // Do not infer vertex count from face count: invalid/non-spherical boundaries
    // must still reach the same nonmanifold and unpaired-edge checks.
    bool dense = p.faces.size() <= 60;
    if (dense) {
        BWScratchTable localIds;
        int vertexCount = 0;
        p.denseVertices.resize(p.faces.size());
        for (size_t i = 0; i < p.faces.size() && dense; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (p.faces[i].vertices[j] < 0) { dense = false; break; }
                int64_t& id = localIds[p.faces[i].vertices[j]];
                if (!id) id = ++vertexCount;
                if (id > 32) { dense = false; break; }
                p.denseVertices[i][j] = static_cast<uint8_t>(id - 1);
            }
        }
    }
    if (dense) {
        // Epoch tags clear only touched entries; 8 bits encode face*3+side+1,
        // bit 8 marks a paired edge. Larger boundaries use the existing hash.
        if (++p.denseEpoch == (1u << 23)) { p.denseEdges.fill(0); p.denseEpoch = 1; }
        const uint32_t epoch = p.denseEpoch << 9;
        size_t unpaired = 0;
        for (size_t i = 0; i < p.faces.size(); ++i) for (int side = 0; side < 3; ++side) {
            int a = p.denseVertices[i][(side + 1) % 3], b = p.denseVertices[i][(side + 2) % 3];
            if (a > b) std::swap(a, b);
            uint32_t& entry = p.denseEdges[a * 32 + b];
            if ((entry & ~511u) != epoch) {
                entry = epoch | static_cast<uint32_t>(i * 3 + side + 1); ++unpaired;
            } else {
                if (entry & 256u) return false;
                const unsigned previous = (entry & 255u) - 1;
                const int face = previous / 3, previousSide = previous % 3;
                entry |= 256u; --unpaired;
                p.faces[i].adjacent[side] = face;
                p.faces[i].adjacentFace[side] = previousSide + 1;
                p.faces[face].adjacent[previousSide] = static_cast<int>(i);
                p.faces[face].adjacentFace[previousSide] = side + 1;
            }
        }
        return unpaired == 0;
    }
    p.edges.clear();
    size_t unpaired = 0;
    for (int i = 0; i < static_cast<int>(p.faces.size()); ++i) {
        const auto vertices = p.faces[i].vertices;
        for (int side = 0; side < 3; ++side) {
            int a = vertices[(side + 1) % 3], b = vertices[(side + 2) % 3];
            if (a > b) std::swap(a, b);
            const uint64_t key = (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
            bool inserted = false;
            BWEdgeTable::Slot& previous = p.edges.insert(key, i, side, inserted);
            if (inserted) ++unpaired;
            else {
                if (previous.paired) return false;
                previous.paired = true;
                --unpaired;
                p.faces[i].adjacent[side] = previous.face;
                p.faces[i].adjacentFace[side] = previous.side + 1;
                p.faces[previous.face].adjacent[previous.side] = i;
                p.faces[previous.face].adjacentFace[previous.side] = side + 1;
            }
        }
    }
    if (unpaired != 0) return false;
    return true;
}

BWStatus prepareBWFill(DT& m, BWPlan& p) {
    if (p.status != BWStatus::Ready) return p.status;
    if (p.request.info == 3) {
        if (p.request.initShell > static_cast<int>(p.workingCavity.size())) return p.status = BWStatus::Rejected;
        for (int i = 0; i < p.request.initShell; ++i)
            if (p.workingCavity[i] == -1) return p.status = BWStatus::Rejected;
    }
    for (size_t i = 0; i < p.workingCavity.size(); ++i) {
        const int tet = p.workingCavity[i];
        if (tet == -1) continue;
        p.cavity.push_back(tet);
        if (p.request.trackAccess) appendWriteTet(p, tet, p.tetInfo[tet]);
        for (int face = 0; face < 4; ++face) if (p.boundaryFaces[i] & (1 << face)) {
            BWFace f;
            f.outside = m.getNeig(tet, face); f.outsideFace = m.getNeigOrd(tet, face);
            if (f.outside >= 0) {
                if (p.request.trackAccess) {
                    int64_t* outsideInfo = checkedTetInfo(m, p, f.outside);
                    if (!outsideInfo) return p.status;
                    appendWriteTet(p, f.outside, *outsideInfo);
                } else if (!validTet(m, f.outside)) return p.status = BWStatus::Stale;
            }
            int ia, ib, ic, id;
            DNC(face, ia, ib, ic, id);
            int b = m.Elems[tet].form[ib], c = m.Elems[tet].form[ic], d = m.Elems[tet].form[id];
            if (b == m.ghost) { std::swap(b, d); std::swap(b, c); }
            else if (c == m.ghost) { std::swap(c, d); std::swap(c, b); }
            f.vertices = {{b, c, d}}; f.geo = m.Elems[tet].geo;
            if (p.request.trackAccess && !p.request.deferP2T)
                for (int v : f.vertices) if (v != m.ghost) appendWriteNode(p, v);
            p.faces.push_back(f);
        }
    }
    if (p.faces.empty()) return p.status = BWStatus::Rejected;
    if (p.request.trackAccess && p.request.node >= 0 && p.request.node != m.ghost)
        appendWriteNode(p, p.request.node);

    if (!pairBWBoundary(p)) return p.status = BWStatus::Rejected;
    p.newElements.resize(p.faces.size(), -1);
    return p.status;
}

BWStatus planBW(DT& m, const BWRequest& request, BWPlan& p) {
    p.request = request;
    if (findBWCavity(m, p) != BWStatus::Ready) return p.status;
    if (adjustBWCavity(m, p) != BWStatus::Ready) return p.status;
    return prepareBWFill(m, p);
}

void commitBW(DT& m, BWPlan& p, int node, const std::vector<int>& slots) {
    if (p.status != BWStatus::Ready || slots.size() != p.faces.size() || p.newElements.size() != slots.size())
        throw std::logic_error("BW commit requires a prepared plan and preallocated slots");
    p.request.node = node;
    std::copy(slots.begin(), slots.end(), p.newElements.begin());
    for (size_t i = 0; i < p.faces.size(); ++i) {
        const BWFace& f = p.faces[i];
        Elem element(node, f.vertices[0], f.vertices[1], f.vertices[2]);
        element.geo = f.geo;
        if (f.outside >= 0) element.neig[0] = (int64_t(f.outside) << 2) | f.outsideFace;
        for (int side = 0; side < 3; ++side)
            element.neig[side + 1] = (int64_t(slots[f.adjacent[side]]) << 2) | f.adjacentFace[side];
        m.Elems[slots[i]] = element;
        if (f.outside >= 0) m.Elems[f.outside].neig[f.outsideFace] = int64_t(slots[i]) << 2;
        m.Nodes[node].tet = slots[i];
        for (int v : f.vertices) {
            if (v == m.ghost) p.ghostTet = slots[i];
            else if (!p.request.deferP2T) m.Nodes[v].tet = slots[i];
        }
    }
    p.status = BWStatus::Inserted;
}

void finishBW(DT& m, BWPlan& p, int thread) {
    if (p.status != BWStatus::Inserted) return;
    if (m.QuantityControl) {
        for (const BWFace& f : p.faces) ++m.GeoNum[f.geo];
        for (int t : p.cavity) --m.GeoNum[m.Elems[t].geo];
    }
    for (int t : p.cavity) {
        if (thread < 0) m.DelEle(t); else m.DelEle(t, thread);
    }
    if (p.request.deferP2T) {
        // All independent commits have finished. Shared old vertices may now
        // choose any incident new tetrahedron without racing another worker.
        for (size_t i = 0; i < p.faces.size(); ++i)
            for (int v : p.faces[i].vertices)
                if (v != m.ghost) m.Nodes[v].tet = p.newElements[i];
    }
    if (p.ghostTet >= 0) m.Nodes[m.ghost].tet = p.ghostTet;
    if (m.improve_step) for (int t : p.newElements) if (!m.isvirtualtet(t)) m.updateQuality(t);
}
}
