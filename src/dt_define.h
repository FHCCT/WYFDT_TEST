#ifndef dt_define_h
#define dt_define_h

#include <atomic>
#include <vector>

namespace dt {
#define PI 3.1415926535897932
#define ANGLE2RADIO(x) (x*PI)/180
#define RADIO2ANGLE(x) (x*180.0)/PI
#define EXCEPTIONSTRING(x) std::logic_error((x).c_str())

	struct Node
	{
		double pt[3];			/* coords. */
		double space;			/* space */
		int tet;				/* point to tet*/
		int info;				/* for B-W point insertion*/
		int type;				/* 0:
								*  1:boudary vertex == Facet vertex
								*  3:Segment vertex
								*  else:Corner vertex
								*/
		std::atomic<int> occupying;

		Node() {
			pt[0] = pt[1] = pt[2] = 0;
			type = tet = space = info = 0; occupying.store(-1);
		};
		Node(double x, double y, double z) {
			pt[0] = x; pt[1] = y; pt[2] = z;
			type = tet = space = info = 0; occupying.store(-1);
		};
		Node(double x, double y, double z, double s) {
			pt[0] = x; pt[1] = y; pt[2] = z; space = s;
			type = tet = info = 0; occupying.store(-1);
		};
		// 自定义拷贝构造函数，正确拷贝 std::atomic
		Node(const Node& other) {
			pt[0] = other.pt[0]; pt[1] = other.pt[1]; pt[2] = other.pt[2];
			space = other.space; tet = other.tet; info = other.info; type = other.type;
			occupying.store(other.occupying.load());  // 通过 load() 进行拷贝
		}

		Node& operator=(const Node& other) {
			if (this == &other) {
				return *this;
			}
			pt[0] = other.pt[0]; pt[1] = other.pt[1]; pt[2] = other.pt[2];
			space = other.space; tet = other.tet; info = other.info; type = other.type;
			occupying.store(other.occupying.load());  // 通过 load() 进行拷贝
			return *this;
		}
	};

	struct Elem
	{
		int form[4];
		int64_t neig[4];		/* neig << 4 + order*/
		int64_t info;			/* bit 0, 1, for point insert,
								 * 0 for test flag in remove outer tet,
								 * 31 for find sphere
								 * 31 for isMeshEdge/Face
								 * 30,29 for point insert
								 * 28 for meshrefine,if it is tried
								 * 0 for finddirection,boundary recover
								 * all for this tet in how much stars */
		int geo;				// 0 is virtual
		double q;				//quality

		Elem() {
			form[0] = form[1] = form[2] = form[3] = -1;
			neig[0] = neig[1] = neig[2] = neig[3] = -1;
			info = 0, geo = -1, q = -1;
		};
		Elem(int a, int b, int c, int d) {
			form[0] = a; form[1] = b; form[2] = c; form[3] = d;
			neig[0] = neig[1] = neig[2] = neig[3] = -1;
			info = 0, geo = -1, q = -1;
		};
		Elem(const Elem& other) {
			form[0] = other.form[0]; form[1] = other.form[1]; form[2] = other.form[2]; form[3] = other.form[3];
			neig[0] = other.neig[0]; neig[1] = other.neig[1]; neig[2] = other.neig[2]; neig[3] = other.neig[3];
			info = other.info, geo = other.geo, q = other.q;
		};
	};

	struct SurTri
	{
		int form[3];			  /* forming points */
		//int edgs[3];			  /* surface edges */
		int parent;				  /*parent */
		int info;				  /* 0=init;
								   * 1=recovered;
								   * < 0 flip deepth;(between boundary recover)
								   * -100 ,destroyed surEdg
								   * else:subTri */
	};

	struct SurEdg
	{
		int iStart, iEnd;		/* start & end point */
		std::vector<int> face;
		//int parent;
		int constrain;			//1,2(given)
		int geo;
		int info;				/* 0=init;
								 * 1=recovered;
								 * < 0 flip deepth;(between boundary recover)
								 * -100 ,destroyed surEdg
								 * else:subEdg */

		SurEdg() : iStart(0), iEnd(0), face(), constrain(0), geo(0),info(0) {}
	};

	/*coding& decoding of elemental entities(nodes / edges / faces)*/
	/*
	* DNC : Decoding Node Codes
	* Known : nc, code of node A
	* Unknown : ib, ic, id, codes of three other nodes
	*/
#define DNC(nc, ia, ib, ic, id)			   \
	ia = (nc);							   \
	switch ((nc))					       \
	{									   \
	case 0:	ib = 1; ic = 2; id = 3;	break; \
	case 1: ib = 3; ic = 2; id = 0; break; \
	case 2:	ib = 0; ic = 1; id = 3;	break; \
	case 3:	ib = 2; ic = 1; id = 0;	break; \
	}
	/*
	* DDNC : Decoding Node Codes
	* Known : ic,id
	* Unknown : ia, ib
	*/
#define DDNC(ia, ib, ic, id)				\
	switch ((ic)){							\
	case 0:									\
		switch((id)){						\
		case 1:ia = 2; ib = 3;	break;		\
		case 2:ia = 3; ib = 1;	break;		\
		case 3:ia = 1; ib = 2;	break;		\
		}break;								\
	case 1:									\
		switch((id)){						\
		case 0:ia = 3; ib = 2;	break;		\
		case 2:ia = 0; ib = 3;	break;		\
		case 3:ia = 2; ib = 0;	break;		\
		}break;								\
	case 2:									\
		switch((id)){						\
		case 0:ia = 1; ib = 3;	break;		\
		case 1:ia = 3; ib = 0;	break;		\
		case 3:ia = 0; ib = 1;	break;		\
		}break;								\
	case 3:									\
		switch((id)){						\
		case 0:ia = 2; ib = 1;	break;		\
		case 1:ia = 0; ib = 2;	break;		\
		case 2:	ia = 1; ib = 0;	break;		\
		}break;								\
	}

	/*
	* DFC: Decoding Face Codes
	* Known : fa, code of face BCD
	* Unknown : ia, ib, ic, id, codes of 4 nodes
	*/
#define DFC(fa, ia, ib, ic, id)				\
	ia = (fa);								\
	switch ((fa))							\
	{										\
	case 0: ib = 1; ic = 3; id = 2; break;	\
	case 1:	ib = 2; ic = 3; id = 0; break;	\
	case 2:	ib = 0; ic = 3; id = 1; break;	\
	case 3:	ib = 0; ic = 1; id = 2; break;	\
	}

	//Edge point ID
	const int Egid[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
	// multiscale_sort Table
	const int  Trans[8][3][8] = {
		{{0,2,6,4,5,7,3,1 },{0,4,5,1,3,7,6,2 },{0,1,3,2,6,7,5,4 }},
		{{1,3,7,5,4,6,2,0 },{1,5,4,0,2,6,7,3 },{1,0,2,3,7,6,4,5 }},
		{{2,0,4,6,7,5,1,3 },{2,6,7,3,1,5,4,0 },{2,3,1,0,4,5,7,6 }},
		{{3,1,5,7,6,4,0,2 },{3,7,6,2,0,4,5,1 },{3,2,0,1,5,4,6,7 }},
		{{4,6,2,0,1,3,7,5 },{4,0,1,5,7,3,2,6 },{4,5,7,6,2,3,1,0 }},
		{{5,7,3,1,0,2,6,4 },{5,1,0,4,6,2,3,7 },{5,4,6,7,3,2,0,1 }},
		{{6,4,0,2,3,1,5,7 },{6,2,3,7,5,1,0,4 },{6,7,5,4,0,1,3,2 }},
		{{7,5,1,3,2,0,4,6 },{7,3,2,6,4,0,1,5 },{7,6,4,5,1,0,2,3 }}
	};
	const int HbTab[8] = { 0,1,0,2,0,1,0,0 };
}
#endif