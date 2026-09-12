#ifndef dt_hash
#define dt_hash

#include <functional>
#include <math.h>
#include <map>
#include <cassert>
#include <utility>
#include <algorithm>
#include <cmath>
#include <vector>

/*************************** Tri Hash ******************************/
template <class T>
class TriHasher
{
public:
	TriHasher();
	void add(int n1, int n2, int n3, const T& val);
	// Does not overwrite existing values. The pointer is invalidated by bucket mutation.
	std::pair<T*, bool> try_emplace(int n1, int n2, int n3, const T& val);
	T get(int n1, int n2, int n3);
	T* find(int n1, int n2, int n3);
	bool erase(int n1, int n2, int n3);
	void clear() {
		tris.clear();
	}
	void clearNode(int n1) {
		tris[n1].clear();
	}

	//用于将索引整合为long long，加速索引相等比较
	typedef union {
		long long hash;
		struct MyStruct
		{
			int n3;
			int n2;
		}idx;
	}TriIndexUnion;

	struct TriInfo
	{
		TriIndexUnion idxunion;
		T val;
	};

	class Iterator {
	public:
		T& operator*() {
			return it->val;
		}

		T* operator->() {
			return &(it->val);
		}

		bool operator==(const Iterator& rhs) {
			return n1 == rhs.n1 && it == rhs.it;
		}

		bool operator!=(const Iterator& rhs) {
			return n1 != rhs.n1 || it != rhs.it;
		}

		int n1, n2, n3;

	private:
		Iterator(int n1, int n2, int n3, const typename std::vector<TriInfo>::iterator& ait)
			:n1(n1), n2(n2), n3(n3), it(ait) {}

		friend TriHasher;
		typename std::vector<TriInfo>::iterator it;
	};

	Iterator begin();
	Iterator end();
	Iterator next(const Iterator& it);

private:
	std::vector<std::vector<TriInfo>> tris;

	//按哈希顺序重排列节点
	void sort(int& n0, int& n1, int& n2) {
		//保证首地址在三个节点号中均匀分布，而不是只取最小节点为首地址
		int idx = (n0 + n1 + n2) % 3;
		switch (idx)
		{
		case 0:
			//从小到大排列
			if (n0 > n1)
				std::swap(n0, n1);
			if (n1 > n2)
				std::swap(n1, n2);
			if (n0 > n1)
				std::swap(n0, n1);
			break;
		case 1:
			//最小在末位，第二大在首位
			if (n0 < n2)
				std::swap(n0, n2);
			if (n1 < n2)
				std::swap(n1, n2);
			if (n0 > n1)
				std::swap(n0, n1);
			break;
		case 2:
			//最小在第二位，第二大在末位
			if (n0 < n1)
				std::swap(n0, n1);
			if (n2 < n1)
				std::swap(n1, n2);
			if (n2 > n0)
				std::swap(n0, n2);
			break;
		}
	}
};

template <class T>
TriHasher<T>::TriHasher()
{
	tris.emplace_back();
}

template<class T>
T TriHasher<T>::get(int n1, int n2, int n3)
{
	sort(n1, n2, n3);
	const auto& tris_n1 = tris[n1];
	for (const auto& tri : tris_n1) {
		if (tri.idxunion.idx.n2 == n2 && tri.idxunion.idx.n3 == n3) {
			return tri.val;
		}
	}
	return -1;
}

template<class T>
T* TriHasher<T>::find(int n1, int n2, int n3) {
	sort(n1, n2, n3);
	if (n1 >= tris.size())
		return nullptr;
	auto& tris_n1 = tris[n1];
	TriIndexUnion idxunion;
	idxunion.idx.n2 = n2;
	idxunion.idx.n3 = n3;
	for (auto& tri : tris_n1) {
		if (tri.idxunion.hash == idxunion.hash)
			return &tri.val;
	}
	return nullptr;
}

template<class T>
std::pair<T*, bool> TriHasher<T>::try_emplace(int n1, int n2, int n3, const T& val)
{
	sort(n1, n2, n3);
	if (n1 < tris.size()) {
		for (auto& entry : tris[n1]) {
			if (entry.idxunion.idx.n2 == n2 && entry.idxunion.idx.n3 == n3)
				return { &entry.val, false };
		}
	}
	else {
		tris.resize(n1 + 1);
	}
	tris[n1].emplace_back();
	tris[n1].back().idxunion.idx.n2 = n2;
	tris[n1].back().idxunion.idx.n3 = n3;
	tris[n1].back().val = val;
	return { &tris[n1].back().val, true };
}

template<class T>
void TriHasher<T>::add(int n1, int n2, int n3, const T& val)
{
	auto result = try_emplace(n1, n2, n3, val);
	if (!result.second)
		*result.first = val;
}

template <class T>
bool TriHasher<T>::erase(int n1, int n2, int n3) {
	sort(n1, n2, n3);
	if (n1 >= tris.size())
		return false;
	auto& tris_n1 = tris[n1];
	TriIndexUnion idxunion;
	idxunion.idx.n2 = n2;
	idxunion.idx.n3 = n3;
	for (auto it = tris_n1.begin(); it != tris_n1.end(); it++) {
		if (it->idxunion.hash == idxunion.hash) {
			tris_n1.erase(it);
			return true;
		}
	}
	return false;
}

template<class T>
inline typename TriHasher<T>::Iterator TriHasher<T>::begin()
{
	int n = 0;
	for (; n < tris.size() && tris[n].empty(); n++);
	if (n == tris.size())
		return end();
	return Iterator(n, tris[n][0].idxunion.idx.n2, tris[n][0].idxunion.idx.n3, tris[n].begin());
}

template <class T>
inline typename TriHasher<T>::Iterator TriHasher<T>::end()
{
    if (tris.empty()) return Iterator(-1, -1, -1, typename std::vector<TriInfo>::iterator{});
	return Iterator(tris.size() - 1, -1, -1, tris.back().end());
}

template<class T>
inline typename TriHasher<T>::Iterator TriHasher<T>::next(const Iterator& it)
{
	int n = it.n1;
	typename std::vector<TriInfo>::iterator vit = it.it;
	if (n == tris.size() - 1 && vit == tris.back().end())
		return it;
	vit++;
	while (true) {
		if (n == tris.size() - 1 && vit == tris.back().end())
			return end();
		if (vit == tris[n].end()) {
			n++;
			vit = tris[n].begin();
		}
		else
			break;
	}
	return Iterator(n, vit->idxunion.idx.n2, vit->idxunion.idx.n3, vit);
}

/*************************** Edge hash ******************************/

template <class T>
class EdgeHasher {
public:
	EdgeHasher();
	void add(int n1, int n2, const T& val);
	// Does not overwrite existing values. The pointer is invalidated by bucket mutation.
	std::pair<T*, bool> try_emplace(int n1, int n2, const T& val);
	T get(int n1, int n2);
	T* find(int n1, int n2);
	bool erase(int n1, int n2);
	void clear() {
		egs.clear();
		egs.emplace_back();
	}
	void clearNode(int n1) {
		egs[n1].clear();
	}

	struct EdgeInfo {
		int n2;
		T val;
	};

	class Iterator {
	public:
		T& operator*() {
			return it->val;
		}

		T* operator->() {
			return &(it->val);
		}

		bool operator==(const Iterator& rhs) {
			return n1 == rhs.n1 && it == rhs.it;
		}

		bool operator!=(const Iterator& rhs) {
			return n1 != rhs.n1 || it != rhs.it;
		}

		int n1, n2;
	private:
		Iterator(int n1, int n2, const typename std::vector<EdgeInfo>::iterator& ait)
			:n1(n1), n2(n2), it(ait) {}

		friend EdgeHasher;
		typename std::vector<EdgeInfo>::iterator it;
	};

	Iterator begin();
	Iterator end();
	Iterator next(const Iterator& it);

private:
	std::vector<std::vector<EdgeInfo>> egs;

	void sort(int& n1, int& n2) {
		if (n1 > n2)
			std::swap(n1, n2);
	}
};

template <class T>
typename EdgeHasher<T>::Iterator EdgeHasher<T>::begin()
{
	int n = 0;
	for (; n < egs.size() && egs[n].empty(); n++);
	if (n == egs.size())
		return end();
	return Iterator(n, egs[n][0].n2, egs[n].begin());
}

template <class T>
typename EdgeHasher<T>::Iterator EdgeHasher<T>::end()
{
	return Iterator(egs.size() - 1, -1, egs.back().end());
}

template <class T>
inline typename EdgeHasher<T>::Iterator EdgeHasher<T>::next(const EdgeHasher<T>::Iterator& it)
{
	int n = it.n1;
	typename std::vector<EdgeInfo>::iterator vit = it.it;
	if (n == egs.size() - 1 && vit == egs.back().end())
		return it;
	vit++;
	while (true) {
		if (n == egs.size() - 1 && vit == egs.back().end())
			return end();
		if (vit == egs[n].end()) {
			n++;
			vit = egs[n].begin();
		}
		else
			break;
	}
	return Iterator(n, vit->n2, vit);
}

template <class T>
bool EdgeHasher<T>::erase(int n1, int n2)
{
	sort(n1, n2);
	if (n1 >= egs.size())
		return false;
	auto& egs_n1 = egs[n1];
	for (auto it = egs_n1.begin(); it != egs_n1.end(); it++) {
		if (it->n2 == n2) {
			egs_n1.erase(it);
			return true;
		}
	}
	return false;
}

template<class T>
std::pair<T*, bool> EdgeHasher<T>::try_emplace(int n1, int n2, const T& val)
{
	sort(n1, n2);
	if (n1 < egs.size()) {
		for (auto& entry : egs[n1]) {
			if (entry.n2 == n2)
				return { &entry.val, false };
		}
	}
	else {
		egs.resize(n1 + 1);
	}
	egs[n1].emplace_back();
	egs[n1].back().n2 = n2;
	egs[n1].back().val = val;
	return { &egs[n1].back().val, true };
}

template<class T>
void EdgeHasher<T>::add(int n1, int n2, const T& val)
{
	auto result = try_emplace(n1, n2, val);
	if (!result.second)
		*result.first = val;
}

template<class T>
EdgeHasher<T>::EdgeHasher() {
	egs.emplace_back();
}

template<class T>
T EdgeHasher<T>::get(int n1, int n2)
{
	sort(n1, n2);
	const auto& egs_n1 = egs[n1];
	for (const auto& eg : egs_n1) {
		if (eg.n2 == n2)
			return eg.val;
	}
	return -1;
}

template<class T>
T* EdgeHasher<T>::find(int n1, int n2) {
	sort(n1, n2);
	if (n1 >= egs.size())
		return nullptr;
	auto& egs_n1 = egs[n1];
	if (egs_n1.empty())
		return nullptr;
	for (auto& eg : egs_n1) {
		if (eg.n2 == n2)
			return &eg.val;
	}
	return nullptr;
}

#endif