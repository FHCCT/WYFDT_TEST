/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 面向网格生成的几何函数库 (版本号：0.1)
 * Geometry Function Library for Mesh Generation (Version 0.1)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2012年11月28日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 28/11/2012
 *
 * 联系方式
 *   电话：+86-571-87951883
 *   传真：+86-571-87953167
 *   邮箱：chenjj@zju.edu.cn
 * For further information, please conctact
 *  Tel: +86-571-87951883
 *  Fax: +86-571-87953167
 * Mail: chenjj@zju.edu.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#include <algorithm>
#include <exception>
#include <stdexcept>
#include <math.h>
#include <assert.h>
#include "geom_func.h"

namespace dt {
	namespace GEOM_FUNC
	{

#define EPS_ZERO_ORIENT3D (1.0e-20)
		enum { DEG, NOD, EDG, FAC };
		/*
		 * Add by Wu Bin
		 */
		 /* Compare two face coplanar-tri was in same direction */
		int cmpTwoTriDir(double*, double*, double*, double*, double*, double*);

		// Labels that signify the result of triangle-triangle intersection test.
		//   Two triangles are DISJOINT, or adjoint at a vertex SHAREVERTEX, or
		//   adjoint at an edge SHAREEDGE, or coincident SHAREFACE or INTERSECT.
		enum interresult { DISJOINT, SHAREVERTEX, SHAREEDGE, SHAREFACE, INTEREDGE, INTERSECT };
		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// Triangle-triangle intersection test                                       //
		//                                                                           //
		// The triangle-triangle intersection test is implemented with exact arithm- //
		// etic. It exactly tells whether or not two triangles in three dimensions   //
		// intersect.  Before implementing this test myself,  I tried two C codes    //
		// (implemented by Thomas Moeller and Philippe Guigue, respectively), which  //
		// are all public available. However both of them failed frequently. Another //
		// unconvenience is both codes only tell whether or not the two triangles    //
		// intersect without distinguishing the cases whether they exactly intersect //
		// in interior or they just share a vertex or share an edge. The two latter  //
		// cases are acceptable and should return not intersection in TetGen.        //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////
		enum interresult edge_vert_col_inter(double*, double*, double*);
		enum interresult edge_edge_cop_inter(double*, double*, double*, double*, double*);
		enum interresult tri_vert_cop_inter(double*, double*, double*, double*, double*, int*);
		enum interresult tri_edge_cop_inter(double*, double*, double*, double*, double*, double*);
		enum interresult tri_edge_inter_tail(double*, double*, double*, double*, double*,
			double, double, int*, double*);
		enum interresult tri_edge_inter(double*, double*, double*, double*, double*);

		// Geometric predicates

		double incircle(double*, double*, double*, double*);
		double insphere_sos(double*, double*, double*, double*, double*, int, int, int, int, int);
		bool iscollinear(double*, double*, double*, double eps);
		bool iscoplanar(double*, double*, double*, double*, double vol6, double eps);
		bool iscospheric(double*, double*, double*, double*, double*, double vol24, double eps);
		// Geometric quantities calculators.
		double distance(double* p1, double* p2);

		int lnFacInt(double ln1[], double ln2[],
			double fac1[], double fac2[], double fac3[],
			double pnt[], int* val,
			int edgFac1, int edgFac2, int edgFac3,
			double* pt1, double* pt2, double* pt3)
		{
			/*
			 * Modified By Wu Bin
			 */
			double a[3], b[3], c[3], p[3], q[3];
			double s1, s2;
			double s[3] = { 1.0, 1.0, 1.0 };
			double n1[3], n2[3], n3[3];

			int i = 0;

			for (i = 0; i != 3; ++i)
			{
				a[i] = fac1[i];
				b[i] = fac2[i];
				c[i] = fac3[i];
				p[i] = ln1[i];
				q[i] = ln2[i];
#if 1
				if (edgFac1 > 0 && pt1 != NULL)//code != NULL)
					n1[i] = pt1[i];//(m_pNodes[code[0]]).pt[i];
				if (edgFac2 > 0 && pt2 != NULL)//code != NULL)
					n2[i] = pt2[i];//(m_pNodes[code[1]]).pt[i];
				if (edgFac3 > 0 && pt3 != NULL)//code != NULL)
					n3[i] = pt3[i];//(m_pNodes[code[2]]).pt[i];
#endif
			}

			if (edgFac1 > 0)
			{
				if (orient3d(n1, p, q, b) == 0.0)
					s[0] = 0.0;
				if (orient3d(n1, p, q, c) == 0.0)
					s[2] = 0.0;
			}
			if (edgFac2 > 0)
			{
				if (orient3d(n2, p, q, a) == 0.0)
					s[0] = 0.0;
				if (orient3d(n2, p, q, c) == 0.0)
					s[1] = 0.0;
			}
			if (edgFac3 > 0)
			{
				if (orient3d(n3, p, q, a) == 0.0)
					s[2] = 0.0;
				if (orient3d(n3, p, q, b) == 0.0)
					s[1] = 0.0;
			}

			// Test the locations of p and q with respect to ABC.
			s1 = orient3d(a, b, c, p);
			s2 = orient3d(a, b, c, q);

#if 0
			if (fabs(s1) < EPS_ZERO_ORIENT3D)
				s1 = 0.0;
			if (fabs(s2) < EPS_ZERO_ORIENT3D)
				s2 = 0.0;
#endif

			bool coplanarFlag = (s1 == 0.0 && s2 == 0.0);
			if (coplanarFlag)	/* coplanar */
			{
				return -1;
			}
			else
			{
				enum interresult abcpq = tri_edge_inter_tail(a, b, c, p, q, s1, s2, val, s);
				if (abcpq == DISJOINT)
				{
					return -1;
				}
				else if (abcpq == INTEREDGE || abcpq == INTERSECT)
				{
					if (s1 == 0.0)
					{
						/* The intersect point is P*/
						pnt[0] = p[0];
						pnt[1] = p[1];
						pnt[2] = p[2];
					}
					else if (s2 == 0.0)
					{
						/* The intersect point is Q*/
						pnt[0] = q[0];
						pnt[1] = q[1];
						pnt[2] = q[2];
					}
					else
					{
						/* The intersect point is inside segment PQ*/

						s1 = orient3dexact(a, b, c, p);
						s2 = orient3dexact(a, b, c, q);
						double tmp1 = s1, tmp2 = s2;
						if (tmp1 < 0.0)
						{
							tmp1 = -tmp1;
						}
						else
						{
							tmp2 = -tmp2;
						}

						pnt[0] = fixedSplitPoint(tmp1, tmp2, p[0], q[0]);
						pnt[1] = fixedSplitPoint(tmp1, tmp2, p[1], q[1]);
						pnt[2] = fixedSplitPoint(tmp1, tmp2, p[2], q[2]);
					}
					if (abcpq == INTEREDGE)
						return EDG;
					else
						return FAC;
				}
				else
				{
					if (*val == 0)
					{
						/* The intersect point is A*/
						pnt[0] = a[0];
						pnt[1] = a[1];
						pnt[2] = a[2];
					}
					else if (*val == 1)
					{
						/* The intersect point is B*/
						pnt[0] = b[0];
						pnt[1] = b[1];
						pnt[2] = b[2];
					}
					else
					{
						/* The intersect point is C*/
						pnt[0] = c[0];
						pnt[1] = c[1];
						pnt[2] = c[2];
					}
					return NOD;
				}
			}

			/* End of modified */
		}

		int lnFacInt2(double ln1[], double ln2[],
			double fac1[], double fac2[], double fac3[],
			double pnt[], int* val,
			int edgFac1, int edgFac2, int edgFac3,
			double* pt1, double* pt2, double* pt3, bool bEpsilon)
		{
			/*
			 * Modified By Wu Bin
			 */
			double a[3], b[3], c[3], p[3], q[3];
			double s[3] = { 1.0, 1.0, 1.0 };
			double n1[3], n2[3], n3[3];
			int i = 0;
			int isInt, intTyp, intCod;
			double linep[2][3], facep[3][3];
			int nRet = -1; /* -1 表示不相交 */

			for (i = 0; i < 3; ++i)
			{
				a[i] = fac1[i];
				b[i] = fac2[i];
				c[i] = fac3[i];
				p[i] = ln1[i];
				q[i] = ln2[i];

				linep[0][i] = ln1[i];
				linep[1][i] = ln2[i];
				facep[0][i] = fac1[i];
				facep[1][i] = fac2[i];
				facep[2][i] = fac3[i];

#if 1
				if (edgFac1 > 0 && pt1 != NULL)//code != NULL)
					n1[i] = pt1[i];//(m_pNodes[code[0]]).pt[i];
				if (edgFac2 > 0 && pt2 != NULL)//code != NULL)
					n2[i] = pt2[i];//(m_pNodes[code[1]]).pt[i];
				if (edgFac3 > 0 && pt3 != NULL)//code != NULL)
					n3[i] = pt3[i];//(m_pNodes[code[2]]).pt[i];
#endif
			}

			if (edgFac1 > 0)
			{
				if (orient3d(n1, p, q, b) == 0.0)
					s[0] = 0.0;
				if (orient3d(n1, p, q, c) == 0.0)
					s[2] = 0.0;
			}
			if (edgFac2 > 0)
			{
				if (orient3d(n2, p, q, a) == 0.0)
					s[0] = 0.0;
				if (orient3d(n2, p, q, c) == 0.0)
					s[1] = 0.0;
			}
			if (edgFac3 > 0)
			{
				if (orient3d(n3, p, q, a) == 0.0)
					s[2] = 0.0;
				if (orient3d(n3, p, q, b) == 0.0)
					s[1] = 0.0;
			}

#if 1
			isInt = lin_tri_intersect3d(linep, facep, &intTyp, &intCod, pnt, bEpsilon);
#else
			isInt = lin_tri_intersect3d_idx(iline, iface, linep, facep, &intTyp, &intCod, pnt, bEpsilon);
#endif

			switch (intTyp)
			{
			case LTI_INTERSECT_NUL:
				nRet = -1;
				break;
			case LTI_INTERSECT_NOD:
				nRet = NOD;
				*val = intCod;
				break;
			case LTI_INTERSECT_EDG:
				nRet = EDG;
				*val = intCod;
				break;
			case LTI_INTERSECT_FAC:
				nRet = FAC;
				break;
			case LTI_INTERSECT_DEG_FACE:
				printf("Warning in lnFacInt2(...). Degenerate triangle.\n");
			default:
				//	assert(intTyp == LTI_INTERSECT_INS);
				nRet = DEG;
				break;
			}

			return nRet;
		}

		int cmpTwoTriDir(double* A, double* B, double* C, double* D, double* E, double* F)
		{
			// pq is coplanar with abc.  Calculate a point which is exactly not
			//   coplanar with a, b, and c.
			double R[3], N[3];
			double ax, ay, az, bx, by, bz;

			ax = A[0] - B[0];
			ay = A[1] - B[1];
			az = A[2] - B[2];
			bx = A[0] - C[0];
			by = A[1] - C[1];
			bz = A[2] - C[2];
			N[0] = ay * bz - by * az;
			N[1] = az * bx - bz * ax;
			N[2] = ax * by - bx * ay;
			// The normal should not be a zero vector (otherwise, abc are collinear).
#ifdef SELF_CHECK
			assert((fabs(N[0]) + fabs(N[1]) + fabs(N[2])) > 0.0);
#endif
			// The reference point R is lifted from A to the normal direction with
			//   a distance d = average edge length of the triangle abc.
			R[0] = N[0] + A[0];
			R[1] = N[1] + A[1];
			R[2] = N[2] + A[2];
			// Becareful the case: if the non-zero component(s) in N is smaller than
			//   the machine epsilon (i.e., 2^(-16) for double), R will exactly equal
			//   to A due to the round-off error.  Do check if it is.
			if (R[0] == A[0] && R[1] == A[1] && R[2] == A[2])
			{
				int i, j;
				for (i = 0; i < 3; i++)
				{
#ifdef SELF_CHECK
					assert(R[i] == A[i]);
#endif
					j = 2;
					do
					{
						if (N[i] > 0.0)
						{
							N[i] += (j * macheps);
						}
						else
						{
							N[i] -= (j * macheps);
						}
						R[i] = N[i] + A[i];
						j *= 2;
					} while (R[i] == A[i]);
				}
			}

			double s1 = orient3d(A, B, C, R), s2 = orient3d(D, E, F, R);
			if (s1 == 0.0 || s2 == 0.0)
			{
				printf("The new test point R was coplanar with to face!!\n");
				return -1;
			}
			else
			{
				if (s1 < 0.0 && s2 < 0.0)
					return 1;
				else if (s1 > 0.0 && s2 > 0.0)
					return 1;
				else
					return 0;
			}

			return -1;
		}

		/*
		 * Add by Wu Bin
		 */
		 ///////////////////////////////////////////////////////////////////////////////
		 //                                                                           //
		 // edge_vert_col_inter()    Test whether an edge (ab) and a collinear vertex //
		 //                          (p) are intersecting or not.                     //
		 //                                                                           //
		 // Possible cases are p is coincident to a (p = a), or to b (p = b), or p is //
		 // inside ab (a < p < b), or outside ab (p < a or p > b). These cases can be //
		 // quickly determined by comparing the corresponding coords of a, b, and p   //
		 // (which are not all equal).                                                //
		 //                                                                           //
		 // The return value indicates one of the three cases: DISJOINT, SHAREVERTEX  //
		 // (p = a or p = b), and INTERSECT (a < p < b).                              //
		 //                                                                           //
		 ///////////////////////////////////////////////////////////////////////////////

		enum interresult edge_vert_col_inter(double* A, double* B,
			double* P)
		{
			int i = 0;
			do {
				if (A[i] < B[i]) {
					if (P[i] < A[i]) {
						return DISJOINT;
					}
					else if (P[i] > A[i]) {
						if (P[i] < B[i]) {
							return INTERSECT;
						}
						else if (P[i] > B[i]) {
							return DISJOINT;
						}
						else {
							// assert(P[i] == B[i]);
							return SHAREVERTEX;
						}
					}
					else {
						// assert(P[i] == A[i]);
						return SHAREVERTEX;
					}
				}
				else if (A[i] > B[i]) {
					if (P[i] < B[i]) {
						return DISJOINT;
					}
					else if (P[i] > B[i]) {
						if (P[i] < A[i]) {
							return INTERSECT;
						}
						else if (P[i] > A[i]) {
							return DISJOINT;
						}
						else {
							// assert(P[i] == A[i]);
							return SHAREVERTEX;
						}
					}
					else {
						// assert(P[i] == B[i]);
						return SHAREVERTEX;
					}
				}
				// i-th coordinates are equal, try i+1-th;
				i++;
			} while (i < 3);
			// Should never be here.
			return DISJOINT;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// edge_edge_cop_inter()    Test whether two coplanar edges (ab, and pq) are //
		//                          intersecting or not.                             //
		//                                                                           //
		// Possible cases are ab and pq are disjointed, or proper intersecting (int- //
		// ersect at a point other than their endpoints), or both collinear and int- //
		// ersecting, or sharing at a common endpoint, or are coincident.            //
		//                                                                           //
		// A reference point R is required, which is exactly not coplanar with these //
		// two edges.  Since the caller knows these two edges are coplanar, it must  //
		// be able to provide (or calculate) such a point.                           //
		//                                                                           //
		// The return value indicates one of the four cases: DISJOINT, SHAREVERTEX,  //
		// SHAREEDGE, and INTERSECT.                                                 //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////
		enum interresult  edge_edge_cop_inter(double* A, double* B,
			double* P, double* Q, double* R)
		{
			double s1, s2, s3, s4;

#ifdef SELF_CHECK
			assert(R != NULL);
#endif
			s1 = orient3d(A, B, R, P);
			s2 = orient3d(A, B, R, Q);
			if (s1 * s2 > 0.0) {
				// Both p and q are at the same side of ab.
				return DISJOINT;
			}
			s3 = orient3d(P, Q, R, A);
			s4 = orient3d(P, Q, R, B);
			if (s3 * s4 > 0.0) {
				// Both a and b are at the same side of pq.
				return DISJOINT;
			}

			// Possible degenerate cases are:
			//   (1) Only one of p and q is collinear with ab;
			//   (2) Both p and q are collinear with ab;
			//   (3) Only one of a and b is collinear with pq.
			enum interresult abp, abq;
			enum interresult pqa, pqb;

			if (s1 == 0.0) {
				// p is collinear with ab.
				abp = edge_vert_col_inter(A, B, P);
				if (abp == INTERSECT) {
					// p is inside ab.
					return INTERSECT;
				}
				if (s2 == 0.0) {
					// q is collinear with ab. Case (2).
					abq = edge_vert_col_inter(A, B, Q);
					if (abq == INTERSECT) {
						// q is inside ab.
						return INTERSECT;
					}
					if (abp == SHAREVERTEX && abq == SHAREVERTEX) {
						// ab and pq are identical.
						return SHAREEDGE;
					}
					pqa = edge_vert_col_inter(P, Q, A);
					if (pqa == INTERSECT) {
						// a is inside pq.
						return INTERSECT;
					}
					pqb = edge_vert_col_inter(P, Q, B);
					if (pqb == INTERSECT) {
						// b is inside pq.
						return INTERSECT;
					}
					if (abp == SHAREVERTEX || abq == SHAREVERTEX) {
						// either p or q is coincident with a or b.
#ifdef SELF_CHECK
		// ONLY one case is possible, otherwise, shoule be SHAREEDGE.
						assert(abp ^ abq);
#endif
						return SHAREVERTEX;
					}
					// The last case. They are disjointed.
#ifdef SELF_CHECK
					assert((abp == DISJOINT) && (abp == abq && abq == pqa && pqa == pqb));
#endif
					return DISJOINT;
				}
				else {
					// p is collinear with ab. Case (1).
#ifdef SELF_CHECK
					assert(abp == SHAREVERTEX || abp == DISJOINT);
#endif
					return abp;
				}
			}
			// p is NOT collinear with ab.
			if (s2 == 0.0) {
				// q is collinear with ab. Case (1).
				abq = edge_vert_col_inter(A, B, Q);
#ifdef SELF_CHECK
				assert(abq == SHAREVERTEX || abq == DISJOINT || abq == INTERSECT);
#endif
				return abq;
			}

			// We have found p and q are not collinear with ab. However, it is still
			//   possible that a or b is collinear with pq (ONLY one of a and b).
			if (s3 == 0.0) {
				// a is collinear with pq. Case (3).
#ifdef SELF_CHECK
				assert(s4 != 0.0);
#endif
				pqa = edge_vert_col_inter(P, Q, A);
#ifdef SELF_CHECK
				// This case should have been detected in above.
				assert(pqa != SHAREVERTEX);
				assert(pqa == INTERSECT || pqa == DISJOINT);
#endif
				return pqa;
			}
			if (s4 == 0.0) {
				// b is collinear with pq. Case (3).
#ifdef SELF_CHECK
				assert(s3 != 0.0);
#endif
				pqb = edge_vert_col_inter(P, Q, B);
#ifdef SELF_CHECK
				// This case should have been detected in above.
				assert(pqb != SHAREVERTEX);
				assert(pqb == INTERSECT || pqb == DISJOINT);
#endif
				return pqb;
			}

			// ab and pq are intersecting properly.
			return INTERSECT;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// Notations                                                                 //
		//                                                                           //
		// Let ABC be the plane passes through a, b, and c;  ABC+ be the halfspace   //
		// including the set of all points x, such that orient3d(a, b, c, x) > 0;    //
		// ABC- be the other halfspace, such that for each point x in ABC-,          //
		// orient3d(a, b, c, x) < 0.  For the set of x which are on ABC, orient3d(a, //
		// b, c, x) = 0.                                                             //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// tri_vert_copl_inter()    Test whether a triangle (abc) and a coplanar     //
		//                          point (p) are intersecting or not.               //
		//                                                                           //
		// Possible cases are p is inside abc, or on an edge of, or coincident with  //
		// a vertex of, or outside abc.                                              //
		//                                                                           //
		// A reference point R is required. R is exactly not coplanar with abc and p.//
		// Since the caller knows they are coplanar, it must be able to provide (or  //
		// calculate) such a point.                                                  //
		//                                                                           //
		// The return value indicates one of the four cases: DISJOINT, SHAREVERTEX,  //
		// and INTERSECT.                                                            //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		enum interresult tri_vert_cop_inter(double* A, double* B,
			double* C, double* P, double* R, int* val)
		{
			double s1, s2, s3;
			int sign;

#ifdef SELF_CHECK
			assert(R != (double*)NULL);
#endif
			// Adjust the orientation of a, b, c and r, so that we can assume that
			//   r is strictly in ABC- (i.e., r is above ABC wrt. right-hand rule).
			s1 = orient3d(A, B, C, R);
#ifdef SELF_CHECK
			assert(s1 != 0.0);
#endif
			sign = s1 < 0.0 ? 1 : -1;

			// Test starts from here.
			s1 = orient3d(A, B, R, P) * sign;
			if (s1 < 0.0) {
				// p is in ABR-.
				return DISJOINT;
			}
			s2 = orient3d(B, C, R, P) * sign;
			if (s2 < 0.0) {
				// p is in BCR-.
				return DISJOINT;
			}
			s3 = orient3d(C, A, R, P) * sign;
			if (s3 < 0.0) {
				// p is in CAR-.
				return DISJOINT;
			}
			if (s1 == 0.0) {
				// p is on ABR.
				if (s2 == 0.0) {
					// p is on BCR.
#ifdef SELF_CHECK
					assert(s3 > 0.0);
#endif
					// p is coincident with b.
					* val = 1;
					return SHAREVERTEX;
				}
				if (s3 == 0.0) {
					// p is on CAR.
					// p is coincident with a.
					*val = 0;
					return SHAREVERTEX;
				}
				// p is on edge ab.
				*val = 0;
				return INTEREDGE;
			}
			// p is in ABR+.
			if (s2 == 0.0) {
				// p is on BCR.
				if (s3 == 0.0) {
					// p is on CAR.
					// p is coincident with c.
					*val = 2;
					return SHAREVERTEX;
				}
				// p is on edge bc.
				*val = 1;
				return INTEREDGE;
			}
			if (s3 == 0.0) {
				// p is on CAR.
				// p is on edge ca.
				*val = 2;
				return INTEREDGE;
			}

			// p is strictly inside abc.
			return INTERSECT;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// tri_edge_cop_inter()    Test whether a triangle (abc) and a coplanar edge //
		//                         (pq) are intersecting or not.                     //
		//                                                                           //
		// A reference point R is required. R is exactly not coplanar with abc and   //
		// pq.  Since the caller knows they are coplanar, it must be able to provide //
		// (or calculate) such a point.                                              //
		//                                                                           //
		// The return value indicates one of the four cases: DISJOINT, SHAREVERTEX,  //
		// SHAREEDGE, and INTERSECT.                                                 //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		enum interresult tri_edge_cop_inter(double* A, double* B,
			double* C, double* P, double* Q, double* R)
		{
			enum interresult abpq, bcpq, capq;
			enum interresult abcp, abcq;
			int p = 0;

			// Test if pq is intersecting one of edges of abc.
			abpq = edge_edge_cop_inter(A, B, P, Q, R);
			if (abpq == INTERSECT || abpq == SHAREEDGE) {
				return abpq;
			}
			bcpq = edge_edge_cop_inter(B, C, P, Q, R);
			if (bcpq == INTERSECT || bcpq == SHAREEDGE) {
				return bcpq;
			}
			capq = edge_edge_cop_inter(C, A, P, Q, R);
			if (capq == INTERSECT || capq == SHAREEDGE) {
				return capq;
			}

			// Test if p and q is inside abc.
			abcp = tri_vert_cop_inter(A, B, C, P, R, &p);
			if (abcp == INTERSECT) {
				return INTERSECT;
			}
			abcq = tri_vert_cop_inter(A, B, C, Q, R, &p);
			if (abcq == INTERSECT) {
				return INTERSECT;
			}

			// Combine the test results of edge intersectings and triangle insides
			//   to detect whether abc and pq are sharing vertex or disjointed.
			if (abpq == SHAREVERTEX) {
				// p or q is coincident with a or b.
#ifdef SELF_CHECK
				assert(abcp ^ abcq);
#endif
				return SHAREVERTEX;
			}
			if (bcpq == SHAREVERTEX) {
				// p or q is coincident with b or c.
#ifdef SELF_CHECK
				assert(abcp ^ abcq);
#endif
				return SHAREVERTEX;
			}
			if (capq == SHAREVERTEX) {
				// p or q is coincident with c or a.
#ifdef SELF_CHECK
				assert(abcp ^ abcq);
#endif
				return SHAREVERTEX;
			}

			// They are disjointed.
			return DISJOINT;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// tri_edge_inter_tail()    Test whether a triangle (abc) and an edge (pq)   //
		//                          are intersecting or not.                         //
		//                                                                           //
		// s1 and s2 are results of pre-performed orientation tests. s1 = orient3d(  //
		// a, b, c, p); s2 = orient3d(a, b, c, q).  To separate this routine from    //
		// tri_edge_inter() can save two orientation tests in tri_tri_inter().       //
		//                                                                           //
		// The return value indicates one of the four cases: DISJOINT, SHAREVERTEX,  //
		// SHAREEDGE, and INTERSECT.                                                 //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		enum interresult tri_edge_inter_tail(double* A, double* B,
			double* C, double* P, double* Q, double s1, double s2, int* val, double* S)
		{
			double s3, s4, s5;
			int sign;

			if (s1 * s2 > 0.0) {
				// p, q are at the same halfspace of ABC, no intersection.
				return DISJOINT;
			}

			if (s1 * s2 < 0.0) {
				// p, q are both not on ABC (and not sharing vertices, edges of abc).
				// Adjust the orientation of a, b, c and p, so that we can assume that
				//   p is strictly in ABC-, and q is strictly in ABC+.
				sign = s1 < 0.0 ? 1 : -1;
				if (S[0] == 0.0)
					s3 = 0.0;
				else
					s3 = orient3d(A, B, P, Q) * sign;
				if (s3 < 0.0) {
					// q is at ABP-.
					return DISJOINT;
				}
				if (S[1] == 0.0)
					s4 = 0.0;
				else
					s4 = orient3d(B, C, P, Q) * sign;
				if (s4 < 0.0) {
					// q is at BCP-.
					return DISJOINT;
				}
				if (S[2] == 0.0)
					s5 = 0.0;
				else
					s5 = orient3d(C, A, P, Q) * sign;
				if (s5 < 0.0) {
					// q is at CAP-.
					return DISJOINT;
				}

				if (fabs(s3) < EPS_ZERO_ORIENT3D)
					s3 = 0.0;
				if (fabs(s4) < EPS_ZERO_ORIENT3D)
					s4 = 0.0;
				if (fabs(s5) < EPS_ZERO_ORIENT3D)
					s5 = 0.0;

				if (s3 == 0.0) {
					// q is on ABP.
					if (s4 == 0.0) {
						// q is on BCP (and q must in CAP+).
						// pq intersects abc at vertex b.
						*val = 1;
						return SHAREVERTEX;
					}
					if (s5 == 0.0) {
						// q is on CAP (and q must in BCP+).
						// pq intersects abc at vertex a.
						*val = 0;
						return SHAREVERTEX;
					}
					// q in both BCP+ and CAP+.
					// pq crosses ab properly.
					*val = 0;
					return INTEREDGE;
				}
				// q is in ABP+;
				if (s4 == 0.0) {
					// q is on BCP.
					if (s5 == 0.0) {
						// q is on CAP.
						// pq intersects abc at vertex c.
						*val = 2;
						return SHAREVERTEX;
					}
					// pq crosses bc properly.
					*val = 1;
					return INTEREDGE;
				}
				// q is in BCP+;
				if (s5 == 0.0) {
					// q is on CAP.
					// pq crosses ca properly.
					*val = 2;
					return INTEREDGE;
				}
				// q is in CAP+;
				// pq crosses abc properly.
				return INTERSECT;
			}

			if (s1 != 0.0 || s2 != 0.0) {
				// Either p or q is coplanar with abc. ONLY one of them is possible.
				if (s1 == 0.0) {
					// p is coplanar with abc, q can be used as reference point.
					return tri_vert_cop_inter(A, B, C, P, Q, val);
				}
				else {
					// q is coplanar with abc, p can be used as reference point.
					return tri_vert_cop_inter(A, B, C, Q, P, val);
				}
			}
			return DISJOINT;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// tri_edge_inter()    Test whether a triangle (abc) and an edge (pq) are    //
		//                     intersecting or not.                                  //
		//                                                                           //
		// The return value indicates one of the four cases: DISJOINT, SHAREVERTEX,  //
		// SHAREEDGE, and INTERSECT.                                                 //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		enum interresult tri_edge_inter(double* A, double* B,
			double* C, double* P, double* Q)
		{
			double s1, s2;
			double s[3] = { 1.0, 1.0, 1.0 };
			int p = 0;

			// Test the locations of p and q with respect to ABC.
			s1 = orient3d(A, B, C, P);
			s2 = orient3d(A, B, C, Q);

			return tri_edge_inter_tail(A, B, C, P, Q, s1, s2, &p, s);
		}

		///////////////////////////////////////////////////
		//                                                                           //
		// insphere_sos()    Insphere test with symbolic perturbation.               //
		//                                                                           //
		// The input points a, b, c, and d should be non-coplanar. They must be ord- //
		// ered so that they have a positive orientation (as defined by orient3d()), //
		// or the sign of the result will be reversed.                               //
		//                                                                           //
		// Return a positive value if the point e lies inside the circumsphere of a, //
		// b, c, and d; a negative value if it lies outside.                         //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		double insphere_sos(double* pa, double* pb, double* pc, double* pd, double* pe,
			int ia, int ib, int ic, int id, int ie)
		{
			double det;

			det = insphere(pa, pb, pc, pd, pe);
			if (det != 0.0) {
				return det;
			}

			// det = 0.0, use symbolic perturbation.
			double* p[5], * tmpp;
			double sign, det_c, det_d;
			int idx[5], perm, tmp;
			int n, i, j;

			p[0] = pa; idx[0] = ia;
			p[1] = pb; idx[1] = ib;
			p[2] = pc; idx[2] = ic;
			p[3] = pd; idx[3] = id;
			p[4] = pe; idx[4] = ie;

			// Bubble sort the points by the increasing order of the indices.
			n = 5;
			perm = 0; // The number of total swaps.
			for (i = 0; i < n - 1; i++) {
				for (j = 0; j < n - 1 - i; j++) {
					if (idx[j + 1] < idx[j]) {  // compare the two neighbors.
						tmp = idx[j];         // swap idx[j] and idx[j + 1]
						idx[j] = idx[j + 1];
						idx[j + 1] = tmp;
						tmpp = p[j];         // swap p[j] and p[j + 1]
						p[j] = p[j + 1];
						p[j + 1] = tmpp;
						perm++;
					}
				}
			}

			sign = (perm % 2 == 0) ? 1.0 : -1.0;
			det_c = orient3d(p[1], p[2], p[3], p[4]); // orient3d(b, c, d, e)
			if (det_c != 0.0) {
				return sign * det_c;
			}
			det_d = orient3d(p[0], p[2], p[3], p[4]); // orient3d(a, c, d, e)
			return -sign * det_d;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// iscollinear()    Check if three points are approximately collinear.       //
		//                                                                           //
		// 'eps' is a relative error tolerance.  The collinearity is determined by   //
		// the value q = cos(theta), where theta is the angle between two vectors    //
		// A->B and A->C.  They're collinear if 1.0 - q <= epspp.                    //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		bool iscollinear(double* A, double* B, double* C, double eps)
		{
			double abx, aby, abz;
			double acx, acy, acz;
			double Lv, Lw, dd;
			double d, q;

			// This implementation has always used a zero absolute-distance scale.
            // Keep that criterion without process-wide writable storage.
            constexpr double longest = 0.0;
			q = longest * eps;
			q *= q;

			abx = A[0] - B[0];
			aby = A[1] - B[1];
			abz = A[2] - B[2];
			acx = A[0] - C[0];
			acy = A[1] - C[1];
			acz = A[2] - C[2];
			Lv = abx * abx + aby * aby + abz * abz;
			// Is AB (nearly) indentical?
			if (Lv < q) return true;
			Lw = acx * acx + acy * acy + acz * acz;
			// Is AC (nearly) indentical?
			if (Lw < q) return true;
			dd = abx * acx + aby * acy + abz * acz;

			d = (dd * dd) / (Lv * Lw);
			if (d > 1.0) d = 1.0; // Rounding.
			q = 1.0 - sqrt(d); // Notice 0 < q < 1.0.

			return q <= eps;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// iscoplanar()    Check if four points are approximately coplanar.          //
		//                                                                           //
		// 'vol6' is six times of the signed volume of the tetrahedron formed by the //
		// four points. 'eps' is the relative error tolerance.  The coplanarity is   //
		// determined by the value: q = fabs(vol6) / L^3,  where L is the average    //
		// edge length of the tet. They're coplanar if q <= eps.                     //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		bool
			iscoplanar(double* k, double* l, double* m, double* n, double vol6, double eps)
		{
			double L, q;
			double x, y, z;

			if (vol6 == 0.0) return true;

			x = k[0] - l[0];
			y = k[1] - l[1];
			z = k[2] - l[2];
			L = sqrt(x * x + y * y + z * z);
			x = l[0] - m[0];
			y = l[1] - m[1];
			z = l[2] - m[2];
			L += sqrt(x * x + y * y + z * z);
			x = m[0] - k[0];
			y = m[1] - k[1];
			z = m[2] - k[2];
			L += sqrt(x * x + y * y + z * z);
			x = k[0] - n[0];
			y = k[1] - n[1];
			z = k[2] - n[2];
			L += sqrt(x * x + y * y + z * z);
			x = l[0] - n[0];
			y = l[1] - n[1];
			z = l[2] - n[2];
			L += sqrt(x * x + y * y + z * z);
			x = m[0] - n[0];
			y = m[1] - n[1];
			z = m[2] - n[2];
			L += sqrt(x * x + y * y + z * z);
#ifdef SELF_CHECK
			assert(L > 0.0);
#endif
			L /= 6.0;
			q = fabs(vol6) / (L * L * L);

			return q <= eps;
		}

		///////////////////////////////////////////////////////////////////////////////
		//                                                                           //
		// iscospheric()    Check if five points are approximately coplanar.         //
		//                                                                           //
		// 'vol24' is the 24 times of the signed volume of the 4-dimensional simplex //
		// formed by the five points. 'eps' is the relative tolerance. The cosphere  //
		// case is determined by the value: q = fabs(vol24) / L^4,  where L is the   //
		// average edge length of the simplex. They're cosphere if q <= eps.         //
		//                                                                           //
		///////////////////////////////////////////////////////////////////////////////

		bool
			iscospheric(double* k, double* l, double* m, double* n, double* o, double vol24, double eps)
		{
			double L, q;

			// A 4D simplex has 10 edges.
			L = distance(k, l);
			L += distance(l, m);
			L += distance(m, k);
			L += distance(k, n);
			L += distance(l, n);
			L += distance(m, n);
			L += distance(k, o);
			L += distance(l, o);
			L += distance(m, o);
			L += distance(n, o);
#ifdef SELF_CHECK
			assert(L > 0.0);
#endif
			L /= 10.0;
			q = fabs(vol24) / (L * L * L * L);

			return q < eps;
		}

		// distance() computs the Euclidean distance between two points.
		double distance(double* p1, double* p2)
		{
			return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
				(p2[1] - p1[1]) * (p2[1] - p1[1]) +
				(p2[2] - p1[2]) * (p2[2] - p1[2]));
		}

		/* -----------------------------------------------------------------------------------
		 * 辅助函数，假设p1/p2/p3/pi(i=1,2...)共面，将p1/p2/p3/pi向使得p1p2p3投影面积最大的坐标平面投影
		 * 如果p1p2p3是退化面片，返回0，否则，返回1
		 * -----------------------------------------------------------------------------------*/
		int proj_coplanr_points(double p1[3], double p2[3], double p3[3], double p4[][3],
			double proj1[2], double proj2[2], double proj3[2], double proj4[][2], int nOth)
		{
			int i1, i2, j;
			double normalf[3], nx, ny, nz;

			/* 求得平面的法向，沿使得平面投影面积最大的方向投影平面和直线 */
			norm_3p(p1, p2, p3, normalf);
			nx = ((normalf[0] < 0.0) ? -normalf[0] : normalf[0]);
			ny = ((normalf[1] < 0.0) ? -normalf[1] : normalf[1]);
			nz = ((normalf[2] < 0.0) ? -normalf[2] : normalf[2]);

			if (nx > 0.0 || ny > 0.0 || nz > 0.0)
			{
				if ((nx > nz) && (nx >= ny))
				{ //往YOZ平面投影
					i1 = 1;
					i2 = 2;
				}
				else if ((ny > nz) && (ny >= nx))
				{//往ZOX平面投影
					i1 = 2;
					i2 = 0;
				}
				else
				{//往XOY平面投影
					i1 = 0;
					i2 = 1;
				}

				proj1[0] = p1[i1];
				proj1[1] = p1[i2];
				proj2[0] = p2[i1];
				proj2[1] = p2[i2];
				proj3[0] = p3[i1];
				proj3[1] = p3[i2];

				for (j = 0; j < nOth; j++)
				{
					proj4[j][0] = p4[j][i1];
					proj4[j][1] = p4[j][i2];
				}
				return 1;
			}
			else
				return 0;
		}

		/* -----------------------------------------------------------------------------------
		 * 辅助函数，假设p1/p2/p3/p4共面，将p1/p2/p3/p4向使得p1p2p3投影面积最大的坐标平面投影
		 * 如果p1p2p3是退化面片，返回0，否则，返回1
		 * -----------------------------------------------------------------------------------*/
		int proj_four_coplanr_points(double p1[3], double p2[3], double p3[3], double p4[3],
			double proj1[2], double proj2[2], double proj3[2], double proj4[2])
		{
			int i1, i2;
			double normalf[3], nx, ny, nz;

			/* 求得平面的法向，沿使得平面投影面积最大的方向投影平面和直线 */
			norm_3p(p1, p2, p3, normalf);
			nx = ((normalf[0] < 0.0) ? -normalf[0] : normalf[0]);
			ny = ((normalf[1] < 0.0) ? -normalf[1] : normalf[1]);
			nz = ((normalf[2] < 0.0) ? -normalf[2] : normalf[2]);

			if (nx > 0.0 || ny > 0.0 || nz > 0.0)
			{
				if ((nx > nz) && (nx >= ny))
				{ //往YOZ平面投影
					i1 = 1;
					i2 = 2;
				}
				else if ((ny > nz) && (ny >= nx))
				{//往ZOX平面投影
					i1 = 2;
					i2 = 0;
				}
				else
				{//往XOY平面投影
					i1 = 0;
					i2 = 1;
				}

				proj1[0] = p1[i1];
				proj1[1] = p1[i2];
				proj2[0] = p2[i1];
				proj2[1] = p2[i2];
				proj3[0] = p3[i1];
				proj3[1] = p3[i2];
				proj4[0] = p4[i1];
				proj4[1] = p4[i2];

				return 1;
			}
			else
				return 0;
		}

		/* -----------------------------------------------------------------------------------
		 * 辅助函数，假设p1/p2/p3/p4共面，判定直线p1/p4是否包含在内角p2p1p3之内
		 * -----------------------------------------------------------------------------------*/
		int line_cross_tri_coplane3d(double p1[3], double p2[3], double p3[3], double p4[3])
		{
			double q1[2], q2[2], q3[2], q4[2];
			int isCross = 1;

			/* 注意，我们不允许退化成一条直线的平面存在，因此，如果出现这种情形，都直接返回1 */
			if (proj_four_coplanr_points(p1, p2, p3, p4, q1, q2, q3, q4) == 1)
			{/* 投影成功 */
				if (orient2d(q1, q2, q3) > 0.0)
					isCross = orient2d(q1, q2, q4) >= 0.0 && orient2d(q1, q3, q4) <= 0.0;
				else
					isCross = orient2d(q1, q3, q4) >= 0.0 && orient2d(q1, q2, q4) <= 0.0;
			}

			return isCross;
		}

		/* -----------------------------------------------------------------------------------
		 * 这是一个判定1个位于面上的3维点和面的关系
		 * 相交类型：
		 * PNT 相交于1个点（intCod返回点的编号0~2，i代表面的第i个顶点）：
		 * EDG 相交于1条边（intCod返回边的编号0~2，i代表(i,(i+1)%3形成的边)
		 * FAC 相交于1个面
		 * intPnt返回交点的值
		 * -----------------------------------------------------------------------------------*/
		int pnt_tri_intersect3d(double poinp[3], double facep[3][3], int* intTyp, int* intCod, double intPnt[3])
		{
			int isInt = 0;
			double q1[2], q2[2], q3[2], q4[2];
			double jOrt, j1Ort, j2Ort, j3Ort;

			/* 注意，我们不允许退化成一条直线的平面存在，因此，如果出现这种情形，都直接返回1 */
			if (proj_four_coplanr_points(facep[0], facep[1], facep[2], poinp, q1, q2, q3, q4) == 1)
			{/* 投影成功 */
				jOrt = orient2d(q1, q2, q3);
				j1Ort = orient2d(q1, q2, q4);
				j2Ort = orient2d(q2, q3, q4);
				j3Ort = orient2d(q3, q1, q4);

				/* 如果j1Ort/j2Ort/j3Ort有2个等于0 (投影成功，则排除3个都等于0的可能)，则交点在点上 */
				if (j1Ort == 0.0 && j2Ort == 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 1;
				}
				else if (j2Ort == 0.0 && j3Ort == 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 2;
				}
				else if (j3Ort == 0.0 && j1Ort == 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 0;
				}

				if (isInt == 0)
				{
					/* 如果j1Ort/j2Ort/j3Ort只有一个等于0，且另外两个与jOrt同号，则交点在边上 */
					if (j1Ort == 0.0 && j2Ort * jOrt > 0.0 && j3Ort * jOrt > 0.0)
					{
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 0;
					}
					if (j2Ort == 0.0 && j3Ort * jOrt > 0.0 && j1Ort * jOrt > 0.0)
					{
						/* 走到这一步，j3Ort和j1Ort的值必然和jOrt同号 */
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 1;
					}
					if (j3Ort == 0.0 && j1Ort * jOrt > 0.0 && j2Ort * jOrt > 0.0)
					{
						/* 走到这一步，j1Ort和j2Ort的值必然和jOrt同号 */
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 2;
					}
				}

				if (isInt == 0)
				{
					/* 如果j1Ort/j2Ort/j3Ort都与jOrt同号，则交点在面内 */
					if (j1Ort * jOrt > 0.0 && j2Ort * jOrt > 0.0 && j3Ort * jOrt > 0.0)
					{
						isInt = 1;
						*intTyp = LTI_INTERSECT_FAC;
					}
				}

				/* 如果相交，保存交点 */
				if (isInt && intPnt)
				{
					intPnt[0] = poinp[0];
					intPnt[1] = poinp[1];
					intPnt[2] = poinp[2];
				}
			}
			return isInt;
		}
		/**
		 * @brief 这个函数只用于判断是否相交
		 * @param linep 点坐标
		 * @param facep 面坐标
		 * @return int 0--不相交 others--相交（包括所有相交）
		 * @note 这个函数对于浮点误差相当敏感，特别是接近0的项，需要小心修改任何一项，如果后续可能，请将浮点替换为float128（boost）
		 * @author 叶鸿飞
		 */
		int lin_tri_intersect3d_check(double linep[2][3], double facep[3][3]) {
			double i1Ort = orient3d(facep[0], facep[1], facep[2], linep[0]);
			double i2Ort = orient3d(facep[0], facep[1], facep[2], linep[1]);
			double linepProj[2][2], facepProj[3][2], intPntProj[2] = { 0 };
			double j1Ort;
			double j2Ort;
			double j3Ort;

			double EPS = 1e-15;
			if (true)
			{
				if (fabs(i1Ort) < EPS)
				{
					//printf("trunction errors.\n");
					i1Ort = 0.0;
				}
				if (fabs(i2Ort) < EPS)
				{
					i2Ort = 0.0;
					//printf("trunction errors.\n");
				}
			}

			/*我们只关心正负，所以能快则快*/

			// if(abs(i1Ort)<EPS&&abs(i2Ort)<EPS){
			// 	/*这种情况，orient3d并不能处理*/
			// 	/*如果投影后不相交，则一定不相交*/
			// 	double f1[3];
			// 	double f2[3];
			// 	double f3[3];
			// 	double l1[2][3];
			// 	double l2[3];
			// 	for(int i=0;i<3;i++){
			// 		f1[i]=facep[0][i];
			// 		f2[i]=facep[1][i];
			// 		f3[i]=facep[2][i];

			// 		l1[0][i]=linep[0][i];
			// 		l1[1][i]=linep[1][i];
			// 	}
			// 			if (proj_coplanr_points(f1, f2, f3, l1,
			// 		facepProj[0], facepProj[1], facepProj[2], linepProj, 2) == 1)
			// 	{
			// 		int intTyp,intCod;
			// 		int isInt = lin_tri_intersect2d(linepProj, facepProj, &intTyp, &intCod, NULL);
			// 		if (isInt == 0)
			// 		{
			// 			return 0;
			// 		}
			// 	}
			// 	else
			// 	{/* 面片退化了，不处理 */
			// 		printf("Error in lin_tri_intersect3d(...). Degenerate triangle.\n");
			// 		printf("If you use multily normal,please rewrite the function\n");
			// 		return 0;
			// 	}
			// 	/*如果运行到了这里，代表至少有一个点有接近0的浮点误差，而且投影到平面也会相交，运行到这里的概率很小*/
			// 	/*如果只有一个点有浮点误差，那么后续的orient3d就不会出现浮点问题,该问题可留给后面的orint3d判断*/
			// 	/*如果两个点都有浮点误差，那么后续的orint3d也有可能会出问题*/
			// }

			/* ---------------------------------------------------------------------
			 * 直线的端点在面的两侧时，必不相交
			 * --------------------------------------------------------------------*/
			if (i1Ort == 0.0 && i2Ort == 0.0)
			{
				/* -------------------------------------------------------------------
				 * 直线和平面共面时，情况很复杂，直线和面的边界可能有2个交点，此时
				 * 暂时只返回离直线首点近的交点
				 * ----------------------------------------------------------------*/
				if (proj_coplanr_points(facep[0], facep[1], facep[2], linep,
					facepProj[0], facepProj[1], facepProj[2], linepProj, 2) == 1)
				{
					int intTyp, intCod;

					return lin_tri_intersect2d(linepProj, facepProj, &intTyp, &intCod, NULL);

					/* 将交点投影到平面上 */
				}
				else
				{/* 面片退化了，不处理 */
		//			assert(false);
					printf("Error in lin_tri_intersect3d(...). Degenerate triangle.\n");
					return 0;
				}
			}
			else if (i1Ort * i2Ort > 0.0)
			{
				return 0;
			}

			if ((i1Ort == 0.0 && i2Ort != 0.0) || (i1Ort != 0.0 && i2Ort == 0.0))
			{
				/* -------------------------------------------------------------------
				 * 直线只有一个端点在面上时，这个点为交点
				 * ----------------------------------------------------------------*/
				int intEndIdx = i1Ort == 0.0 ? 0 : 1;
				int intTyp, intCod;
				double intPnt[3];
				int isInt = pnt_tri_intersect3d(linep[intEndIdx], facep, &intTyp, &intCod, intPnt);
			}
			else
			{
				/* -------------------------------------------------------------------
				 * 直线2个点都不在平面时，假设(p1,p2,p3)法向指向p4(直线的1个端点，另一个
				 * 负法向端点记为p5，则p4/p5和p1 p2 p3 相交的充要条件是：
				 * p1p2p5p4/p2p3p5p4/p3p1p5p4的体积都大于等于0 (相应orient3d则小于0)
				 * ----------------------------------------------------------------*/
				int p4Idx = i1Ort > 0.0 ? 1 : 0;
				int p5Idx = i1Ort > 0.0 ? 0 : 1;

				j1Ort = orient3d(facep[0], facep[1], linep[p4Idx], linep[p5Idx]);
				j2Ort = orient3d(facep[1], facep[2], linep[p4Idx], linep[p5Idx]);
				j3Ort = orient3d(facep[2], facep[0], linep[p4Idx], linep[p5Idx]);

				/* 如果j1Ort/j2Ort/j3Ort有2个等于0 (因为2个点都不在面上，且面不退化，可排除3个都等于0的可能)，则交点在点上 */
				if (j1Ort == 0.0 && j2Ort == 0.0)
				{
					return 1;
				}
				else if (j2Ort == 0.0 && j3Ort == 0.0)
				{
					return 1;
				}
				else if (j3Ort == 0.0 && j1Ort == 0.0)
				{
					return 1;
				}

				/* 如果j1Ort/j2Ort/j3Ort只有一个等于0，且另外两个都大于0，则交点在边上 */
				if (j1Ort == 0.0 && j2Ort > 0.0 && j3Ort > 0.0)
				{
					return 1;
				}
				else if (j2Ort == 0.0 && j3Ort > 0.0 && j1Ort > 0.0)
				{
					return 1;
				}
				else if (j3Ort == 0.0 && j1Ort > 0.0 && j2Ort > 0.0)
				{
					return 1;
				}
				else if (j1Ort > 0.0 && j2Ort > 0.0 && j3Ort > 0.0)
				{
					return 1;
				}
			}
			return 0;
		}
		/* -----------------------------------------------------------------------------------
		 * 这是一个判定1条线段和三角形是否相交的代码(2D)，并返回交点位置
		 * 相交类型：
		 * PNT 相交于1个点（intCod返回点的编号0~2，i代表面的第i个顶点）：
		 * EDG 相交于1条边（intCod返回边的编号0~2，i代表(i,(i+1)%3形成的边)
		 * FAC 相交于1个面
		 * intPnt返回交点的值
		 * 注意：在共面情形下，一条线可能会和一个面的两条边都相交，此时，
		 * 根据调用该算法的边界恢复过程需求，我们只返回离得最近的那个交点（如果有一个
		 * 交点为直线的顶点，则返回不是顶点的交点）
		 * 当将这个代码用于表面可能相交的情形时，我们需要根据具体情况进行更新
		 * -----------------------------------------------------------------------------------*/
		int lin_tri_intersect2d(double linep[2][2], double facep[3][2], int* intTyp, int* intCod, double intPnt[2])
		{
			int isInt = 0, isInner1, isInner2;
			double iOrt, iOrt1, iOrt2, iOrt3, iOrt11, iOrt12, iOrt13, iOrt21, iOrt22, iOrt23, kOrt1, kOrt2;
			double iOrtArr[3], jOrtArr[3], mOrtArr[3]; /* mOrtArr:面的3个顶点相对于直线的绕向 */
			double w1 = 0.0, w2 = 0.0;
			int m, nodCod, edgCod, candEdgNum = 0, instNodNum = 0, ii;
			double candEdgCod[3], twoIntPnt[2][2], twoIntCod[2], twoIntTyp[2]; /* 可能有2个交点 */
			double d1 = 0.0, d2 = 0.0;
			int iCord;

			/* -------------------------------------------------------------------
			 * 判定2个顶点是否在三角形内部
			 * ----------------------------------------------------------------*/
			iOrt = orient2d(facep[0], facep[1], facep[2]);
			iOrt11 = orient2d(facep[0], facep[1], linep[0]);
			iOrt12 = orient2d(facep[1], facep[2], linep[0]);
			iOrt13 = orient2d(facep[2], facep[0], linep[0]);
			iOrt21 = orient2d(facep[0], facep[1], linep[1]);
			iOrt22 = orient2d(facep[1], facep[2], linep[1]);
			iOrt23 = orient2d(facep[2], facep[0], linep[1]);

			isInner1 = iOrt * iOrt11 >= 0.0 && iOrt * iOrt12 >= 0.0 && iOrt * iOrt13 >= 0.0;
			isInner2 = iOrt * iOrt21 >= 0.0 && iOrt * iOrt22 >= 0.0 && iOrt * iOrt23 >= 0.0;

			*intTyp = LTI_INTERSECT_NUL; /* 预设为不相交 */

			if (isInner1 && isInner2)
			{
				/* 直线两个点都在平面的内部或边界上时，这条直线被包含在平面内部 */
				isInt = 1;
				*intTyp = LTI_INTERSECT_INS;
			}
			else if (isInner1 || isInner2)
			{/* 直线1个点在平面的内部或边界上时，这条直线和平面必有一个交点 */
				iOrtArr[0] = iOrt1 = isInner1 ? iOrt11 : iOrt21;
				iOrtArr[1] = iOrt2 = isInner1 ? iOrt12 : iOrt22;
				iOrtArr[2] = iOrt3 = isInner1 ? iOrt13 : iOrt23;
				jOrtArr[0] = isInner2 ? iOrt11 : iOrt21;
				jOrtArr[1] = isInner2 ? iOrt12 : iOrt22;
				jOrtArr[2] = isInner2 ? iOrt13 : iOrt23;

				/* 如果iOrt1/iOrt2/iOrt3有2个等于0 (投影成功，则排除3个都等于0的可能)，则交点在点上 */
				if (iOrt1 == 0.0 && iOrt2 == 0.0)
				{
					assert(iOrt3 != 0.0);
					*intTyp = LTI_INTERSECT_NOD;
					isInt = 1;
					*intCod = 1;
				}
				else if (iOrt2 == 0.0 && iOrt3 == 0.0)
				{
					assert(iOrt1 != 0.0);
					*intTyp = LTI_INTERSECT_NOD;
					isInt = 1;
					*intCod = 2;
				}
				else if (iOrt3 == 0.0 && iOrt1 == 0.0)
				{
					assert(iOrt2 != 0.0);
					*intTyp = LTI_INTERSECT_NOD;
					isInt = 1;
					*intCod = 0;
				}

				if (isInt == 0)
				{
					/* 如果iOrt1/iOrt2/iOrt3只有一个等于0，且另外两个与iOrt同号，则交点在边上 */
					if (iOrt1 == 0.0)
					{
						/* 走到这一步，iOrt2和iOrt3值必然和iOrt同号 */
						assert(iOrt2 * iOrt > 0.0 && iOrt3 * iOrt > 0.0);
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 0;
					}
					if (iOrt2 == 0.0)
					{
						/* 走到这一步，iOrt3和iOrt1值必然和iOrt同号 */
						assert(iOrt3 * iOrt > 0.0 && iOrt1 * iOrt > 0.0);
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 1;
					}
					if (iOrt3 == 0.0)
					{
						/* 走到这一步，iOrt1和iOrt2值必然和iOrt同号 */
						assert(iOrt1 * iOrt > 0.0 && iOrt2 * iOrt > 0.0);
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 2;
					}
				}

				if (isInt == 1)
				{
					ii = isInner1 ? 0 : 1;
					if (intPnt)
					{
						intPnt[0] = linep[ii][0];
						intPnt[1] = linep[ii][1];
					}
				}

				if (isInt == 0)
				{/* 此时，只有一个交点 */
					mOrtArr[0] = orient2d(linep[0], linep[1], facep[0]);
					mOrtArr[1] = orient2d(linep[0], linep[1], facep[1]);
					mOrtArr[2] = orient2d(linep[0], linep[1], facep[2]);

					/* 第一类情形，面有1个端点在直线上 （走到这，面不可能有2个端点在直线上）*/
					for (nodCod = 0; nodCod < 3; nodCod++)
					{
						if (mOrtArr[nodCod] == 0.0)
						{
							isInt = 1;

							*intTyp = LTI_INTERSECT_NOD;
							*intCod = nodCod;
							if (intPnt)
							{
								intPnt[0] = facep[nodCod][0];
								intPnt[1] = facep[nodCod][1];
							}

							break;
						}
					}

					/* 第二类情形，面有1条边和它相交 */
					if (isInt == 0)
					{
						for (edgCod = 0; edgCod < 3; edgCod++)
						{
							if (mOrtArr[edgCod] * mOrtArr[(edgCod + 1) % 3] < 0.0 && iOrtArr[edgCod] * jOrtArr[edgCod] < 0.0)
							{/* 相交 */
								isInt = 1;

								*intTyp = LTI_INTERSECT_EDG;
								*intCod = edgCod;

								if (intPnt)
								{
									w1 = isInner1 ? jOrtArr[edgCod] : iOrtArr[edgCod];
									w2 = isInner1 ? iOrtArr[edgCod] : jOrtArr[edgCod];
									w1 = fabs(w1) / (fabs(w1) + fabs(w2));
									w2 = 1.0 - w1;
									intPnt[0] = w1 * linep[0][0] + w2 * linep[1][0];
									intPnt[1] = w1 * linep[0][1] + w2 * linep[1][1];
								}
								break;
							}
						}
					}
				}
				else
				{
					isInt = 1;
					/* 显然，直线可能还会和平面有第2个交点 */
					if (*intTyp == LTI_INTERSECT_NOD)
					{/* 判定直线是否在顶点所在的内角上 */
						/* *intCod = 0时，需判定jOrtArr[2]/jOrtArr[0];
						 * *intCod = 1时，需判定jOrtArr[0]/jOrtArr[1];
						 * *intCod = 2时，需判定jOrtArr[1]/jOrtArr[2];
						 */
						assert(*intCod >= 0 && *intCod <= 2);
						nodCod = *intCod;
						if (jOrtArr[nodCod] * iOrt >= 0.0 && jOrtArr[(nodCod + 2) % 3] * iOrt >= 0.0)
						{/* 相交 */
							if (jOrtArr[nodCod] == 0.0)
							{/* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
								/* 边的编号是edgCod, 两个顶点是edgCod和(edgCod+1)%3 */
								edgCod = nodCod;
								*intTyp = LTI_INTERSECT_NOD;
								*intCod = edgCod == nodCod ? (edgCod + 1) % 3 : edgCod;
								if (intPnt)
								{
									intPnt[0] = facep[*intCod][0];
									intPnt[1] = facep[*intCod][1];
								}
							}
							else if (jOrtArr[(nodCod + 2) % 3] == 0.0)
							{/* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
								/* 边的编号是edgCod, 两个顶点是edgCod和(edgCod+1)%3 */
								edgCod = (nodCod + 2) % 3;
								*intTyp = LTI_INTERSECT_NOD;
								*intCod = edgCod == nodCod ? (edgCod + 1) % 3 : edgCod;
								if (intPnt)
								{
									intPnt[0] = facep[*intCod][0];
									intPnt[1] = facep[*intCod][1];
								}
							}
							else
							{/* 交点在边上 */
								/* 当第1个点在平面内部，对应的jOrt代表不在平面点的和第3条边的面积，它越大，
								 * 则直线第1个点的贡献越大
								 */
								edgCod = (nodCod + 1) % 3;
								/* 我们这是只返回第2个顶点的信息，此处需要根据具体应用进行更新 */
								*intTyp = LTI_INTERSECT_EDG;
								*intCod = edgCod;

								if (intPnt)
								{
									w1 = isInner1 ? jOrtArr[edgCod] : iOrtArr[edgCod];
									w2 = isInner1 ? iOrtArr[edgCod] : jOrtArr[edgCod];
									w1 = fabs(w1) / (fabs(w1) + fabs(w2));
									w2 = 1.0 - w1;

									intPnt[0] = w1 * linep[0][0] + w2 * linep[1][0];
									intPnt[1] = w1 * linep[0][1] + w2 * linep[1][1];
								}
							}
						}
					}
					else if (*intTyp == LTI_INTERSECT_EDG)
					{/* 可能和两位两条边相交 */
						assert(*intCod >= 0 && *intCod <= 2);
						edgCod = *intCod;
						/* 另外1个顶点必然也在边的同一侧 */
						if (jOrtArr[edgCod] * iOrt > 0.0)
						{/* 相交 */
							/* 如果另外2条边对应的jOrt值都小于0，还有进一步区分 */
							if (jOrtArr[(edgCod + 1) % 3] * iOrt < 0.0 && jOrtArr[(edgCod + 2) % 3] * iOrt < 0.0)
							{
								/* 判定edgCod + 1的两个端点是否在直线的两侧 */
								edgCod = (edgCod + 1) % 3;
								kOrt1 = orient2d(linep[0], linep[1], facep[edgCod]);
								kOrt2 = orient2d(linep[0], linep[1], facep[(edgCod + 1) % 3]);
								if (kOrt1 == 0.0 || kOrt2 == 0.0)
								{/* 相交于面片的一个顶点 */
									isInt = 1;
									*intTyp = LTI_INTERSECT_NOD;
									*intCod = kOrt1 == 0.0 ? edgCod : (edgCod + 1) % 3;
									if (intPnt)
									{
										intPnt[0] = facep[*intCod][0];
										intPnt[1] = facep[*intCod][1];
									}
									edgCod = -1;
								}
								else if (kOrt1 * kOrt2 >= 0.0)
									edgCod = (edgCod + 2) % 3;
							}
							else if (jOrtArr[(edgCod + 1) % 3] * iOrt < 0.0)
							{
								edgCod = (edgCod + 1) % 3;
							}
							else
							{
								/* 走到这，jOrtArr[(edgCod+2)%3]*iOrt < 0.0 */
								assert(jOrtArr[(edgCod + 2) % 3] * iOrt < 0.0);
								edgCod = (edgCod + 2) % 3;
							}

							if (edgCod >= 0)
							{
								/* 和平面第edgCod条边相交 */
								/* 当第1个点在平面内部，对应的jOrt代表不在平面点的和第3条边的面积，它越大，
								 * 则直线第1个点的贡献越大
								 */

								 /* 我们这是只返回第2个顶点的信息，此处需要根据具体应用进行更新 */
								*intTyp = LTI_INTERSECT_EDG;
								*intCod = edgCod;
								if (intPnt)
								{
									w1 = isInner1 ? jOrtArr[edgCod] : iOrtArr[edgCod];
									w2 = isInner1 ? iOrtArr[edgCod] : jOrtArr[edgCod];
									w1 = fabs(w1) / (fabs(w1) + fabs(w2));
									w2 = 1.0 - w1;

									intPnt[0] = w1 * linep[0][0] + w2 * linep[1][0];
									intPnt[1] = w1 * linep[0][1] + w2 * linep[1][1];
								}
							}
						}
						else if (jOrtArr[edgCod] == 0.0)
						{ /* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
							*intTyp = LTI_INTERSECT_NOD;
							/* 看下一条边的面积值，如果小于0，则是边edgCod的末点 */
							*intCod = jOrtArr[(edgCod + 1) % 3] < 0.0 ? (edgCod + 1) % 3 : edgCod;

							if (intPnt)
							{
								intPnt[0] = facep[*intCod][0];
								intPnt[1] = facep[*intCod][1];
							}
						}
					}
				}/* else (isInt != 0) */
			}/* else if (isInner1 || isInner2) */
			else
			{/* 直线两个点都在外部，仍然可能和平面相交于1个或2个点 */
				iOrtArr[0] = iOrt11;
				iOrtArr[1] = iOrt12;
				iOrtArr[2] = iOrt13;
				jOrtArr[0] = iOrt21;
				jOrtArr[1] = iOrt22;
				jOrtArr[2] = iOrt23;

				if (iOrtArr[0] * jOrtArr[0] <= 0.0 || iOrtArr[1] * jOrtArr[1] <= 0.0 || iOrtArr[2] * jOrtArr[2] <= 0.0)
				{/* 有可能相交 */
					/* 第一类情形，面的某条边是原来直线的子线段 */
					for (edgCod = 0; edgCod < 3; edgCod++)
					{
						if (iOrtArr[edgCod] == 0.0 && jOrtArr[edgCod] == 0.0)
						{
							/* 确定面片这条边的2个点都在边的内部 */
							iCord = fabs(linep[0][0] - linep[1][0]) > fabs(linep[0][1] - linep[1][1]) ? 0 : 1;
							d1 = linep[0][iCord] > linep[1][iCord] ? linep[0][iCord] : linep[1][iCord];
							d2 = linep[0][iCord] > linep[1][iCord] ? linep[1][iCord] : linep[0][iCord];

							if (facep[edgCod][iCord] >= d2 && facep[edgCod][iCord] <= d1)
							{/* 不在延长线上 */
								/* 另外一个点也在线的内部 */
								assert(facep[(edgCod + 1) % 3][iCord] >= d2 && facep[(edgCod + 1) % 3][iCord] <= d1);
								isInt = 1;
								*intTyp = LTI_INTERSECT_EDG;
								*intCod = edgCod;
								/* 只记录远端的那个顶点 */
								if (intPnt)
								{
									if (jOrtArr[(edgCod + 1) % 3] * iOrt < 0.0)
									{/* 记录边的末点，它离直线的起点更远 */
										intPnt[0] = facep[(edgCod + 1) % 3][0];
										intPnt[1] = facep[(edgCod + 1) % 3][1];
									}
									else
									{/* 记录边的首点，它离直线的起点更远 */
										intPnt[0] = facep[edgCod][0];
										intPnt[1] = facep[edgCod][1];
									}
								}

								return 1;
							}
							else
							{
								*intTyp = LTI_INTERSECT_NUL;
								return 0; /* 不相交 */
							}
						}
					}

#if 0
					/* 第二类情形，面的某个端点是直线上一点 */
					if (isInt == 0)
					{
						for (nodCod = 0; nodCod < 3; nodCod++)
						{
							mOrtArr[nodCod] = orient2d(linep[0], linep[1], facep[nodCod]);
							if (mOrtArr[nodCod] == 0.0)
							{
								/* 是否在延长线上 */
								iCord = fabs(linep[0][0] - linep[1][0]) > fabs(linep[0][1] - linep[1][1]) ? 0 : 1;
								d1 = linep[0][iCord] > linep[1][iCord] ? linep[0][iCord] : linep[1][iCord];
								d2 = linep[0][iCord] > linep[1][iCord] ? linep[1][iCord] : linep[0][iCord];

								if (facep[nodCod][iCord] >= d2 && facep[nodCod][iCord] <= d1)
								{/* 不在延长线上 */
									isInt = 1;
									*intTyp = LTI_INTERSECT_NOD;
									*intCod = nodCod;

									intPnt[0] = facep[nodCod][0];
									intPnt[1] = facep[nodCod][1];
									break;
								}
							}
						}
					}
#endif

					/* 第2类情形，面有1个端点在直线上 （走到这，面不可能有2个端点在直线上）*/
					mOrtArr[0] = orient2d(linep[0], linep[1], facep[0]);
					mOrtArr[1] = orient2d(linep[0], linep[1], facep[1]);
					mOrtArr[2] = orient2d(linep[0], linep[1], facep[2]);
					for (nodCod = 0; nodCod < 3; nodCod++)
					{
						if (mOrtArr[nodCod] == 0.0)
						{
							/* 是否在延长线上 */
							iCord = fabs(linep[0][0] - linep[1][0]) > fabs(linep[0][1] - linep[1][1]) ? 0 : 1;
							d1 = linep[0][iCord] > linep[1][iCord] ? linep[0][iCord] : linep[1][iCord];
							d2 = linep[0][iCord] > linep[1][iCord] ? linep[1][iCord] : linep[0][iCord];

							if (facep[nodCod][iCord] < d2 || facep[nodCod][iCord] > d1)
							{/* 延长线上，不相交 */
								*intTyp = LTI_INTERSECT_NUL;
								return 0;
							}

							/* 不在延长线上，第1个交点即为该顶点 */
							isInt = 1;

							instNodNum = 0;
							twoIntTyp[instNodNum] = LTI_INTERSECT_NOD;
							twoIntCod[instNodNum] = nodCod;
							twoIntPnt[instNodNum][0] = facep[nodCod][0];
							twoIntPnt[instNodNum][1] = facep[nodCod][1];
							instNodNum++;

							/* 看看和nodCod相对的边是否和直线相交 */
							edgCod = (nodCod + 1) % 3;
							if (mOrtArr[edgCod] * mOrtArr[(edgCod + 1) % 3] < 0.0 && iOrtArr[edgCod] * jOrtArr[edgCod] < 0.0)
							{
								twoIntTyp[instNodNum] = LTI_INTERSECT_EDG;
								twoIntCod[instNodNum] = edgCod;

								w1 = fabs(jOrtArr[edgCod]);
								w2 = fabs(iOrtArr[edgCod]);
								w1 = w1 / (w1 + w2);
								w2 = 1.0 - w1;
								twoIntPnt[instNodNum][0] = w1 * linep[0][0] + w2 * linep[1][0];
								twoIntPnt[instNodNum][1] = w1 * linep[0][1] + w2 * linep[1][1];
								instNodNum++;
							}

							break;
						}
					}

					/* 第3类情形，面的2条边和直线有交点 */
					if (isInt == 0)
					{
						/* 如果直线的2个端点在面片3条边的某条边不同侧时 */
						for (candEdgNum = 0, edgCod = 0; edgCod < 3; edgCod++)
						{
							/* 注意，此时没有等于0的判断，因为如果等于0，说明一个顶点和这条曲线共线
								* 这和这条边的交点必然也位于在这条边上，这在前面已有判断 */
							if (iOrtArr[edgCod] * jOrtArr[edgCod] < 0.0)
							{
								candEdgCod[candEdgNum++] = edgCod;
							}
						}

						if (candEdgNum > 0)
						{
							for (m = 0, instNodNum = 0; m < candEdgNum; m++)
							{
								edgCod = candEdgCod[m];
								/* 边的两个端点是否异号 */
								if (mOrtArr[edgCod] * mOrtArr[(edgCod + 1) % 3] < 0.0)
								{/* 相交 */
									isInt = 1;
									assert(instNodNum < 2);
									if (instNodNum >= 2)
									{
										//							printf("Error in lin_tri_intersect2d(...). Degenerate triangle.\n");
										*intTyp = LTI_INTERSECT_DEG_FACE;
										return 1; /* 返回相交 */
									}
									twoIntTyp[instNodNum] = LTI_INTERSECT_EDG;
									twoIntCod[instNodNum] = edgCod;

									w1 = fabs(jOrtArr[edgCod]);
									w2 = fabs(iOrtArr[edgCod]);
									w1 = w1 / (w1 + w2);
									w2 = 1.0 - w1;
									twoIntPnt[instNodNum][0] = w1 * linep[0][0] + w2 * linep[1][0];
									twoIntPnt[instNodNum][1] = w1 * linep[0][1] + w2 * linep[1][1];
									instNodNum++;
								}
							}
						}
					}

					if (isInt != 0)
					{
						assert(instNodNum == 1 || instNodNum == 2);
						ii = 0;
						if (instNodNum == 2)
						{/* 取离直线起点最远的交点 */
							d1 = (linep[0][0] - twoIntPnt[0][0]) * (linep[0][0] - twoIntPnt[0][0]) +
								(linep[0][1] - twoIntPnt[0][1]) * (linep[0][1] - twoIntPnt[0][1]);
							d2 = (linep[0][0] - twoIntPnt[1][0]) * (linep[0][0] - twoIntPnt[1][0]) +
								(linep[0][1] - twoIntPnt[1][1]) * (linep[0][1] - twoIntPnt[1][1]);
							ii = d1 < d2 ? 1 : 0;
						}

						*intTyp = twoIntTyp[ii];
						*intCod = twoIntCod[ii];
						if (intPnt)
						{
							intPnt[0] = twoIntPnt[ii][0];
							intPnt[1] = twoIntPnt[ii][1];
						}
					}
				}
			}

			return isInt;
		}

		/* -----------------------------------------------------------------------------------
		 * 这是一个判定1条线段和三角形是否相交的代码(3D)，并返回交点位置
		 * 相交类型：
		 * PNT 相交于1个点（intCod返回点的编号0~2，i代表面的第i个顶点）：
		 * EDG 相交于1条边（intCod返回边的编号0~2，i代表(i,(i+1)%3形成的边)
		 * FAC 相交于1个面
		 * intPnt返回交点的值
		 * -----------------------------------------------------------------------------------*/
		int lin_tri_intersect3d(double linep[2][3], double facep[3][3], int* intTyp, int* intCod, double intPnt[3], bool bEpsilon)
		{
			double i1Ort = 0, i2Ort = 0;
			int isInt = 0;
			int intEndIdx = 0, p4Idx = 0, p5Idx = 0;
			double j1Ort = 0, j2Ort = 0, j3Ort = 0;
			double d1 = 0, d = 0, w1 = 0, w2 = 0;
			double linepProj[2][2] = { 0 }, facepProj[3][2] = {0}, intPntProj[2] = { 0 };

			*intTyp = LTI_INTERSECT_NUL;
			*intCod = -1;

			i1Ort = orient3d(facep[0], facep[1], facep[2], linep[0]);
			i2Ort = orient3d(facep[0], facep[1], facep[2], linep[1]);

			//#WARNING yhf: 这个选项不要开，复杂情况边界恢复时若开了这个选项很容易出现误判
			bEpsilon = false;

			if (bEpsilon)
			{
				if (fabs(i1Ort) < 1.0e-18)
				{
					//printf("trunction errors.\n");
					i1Ort = 0.0;
				}
				if (fabs(i2Ort) < 1.0e-18)
				{
					i2Ort = 0.0;
					//printf("trunction errors.\n");
				}
			}

			/* ---------------------------------------------------------------------
			 * 直线的端点在面的两侧时，必不相交
			 * --------------------------------------------------------------------*/
			if (i1Ort * i2Ort > 0.0)
			{
				*intTyp = LTI_INTERSECT_NUL;
				return 0;
			}

			if (i1Ort == 0.0 && i2Ort == 0.0)
			{
				/* -------------------------------------------------------------------
				 * 直线和平面共面时，情况很复杂，直线和面的边界可能有2个交点，此时
				 * 暂时只返回离直线首点近的交点
				 * ----------------------------------------------------------------*/
				if (proj_coplanr_points(facep[0], facep[1], facep[2], linep,
					facepProj[0], facepProj[1], facepProj[2], linepProj, 2) == 1)
				{
					isInt = lin_tri_intersect2d(linepProj, facepProj, intTyp, intCod, intPnt ? intPntProj : NULL);

					/* 将交点投影到平面上 */
					if (isInt == 1)
					{
						if (*intTyp == LTI_INTERSECT_NOD)
						{
							if (intPnt)
							{
								intPnt[0] = facep[*intCod][0];
								intPnt[1] = facep[*intCod][1];
								intPnt[2] = facep[*intCod][2];
							}
						}
						else if (isInt == 1)
						{
							if (intPnt)
							{
								/* 当我们需要更精确的坐标时，我们需要更小心处理浮点误差 */

								d1 = (linepProj[0][0] - intPntProj[0]) * (linepProj[0][0] - intPntProj[0]) +
									(linepProj[0][1] - intPntProj[1]) * (linepProj[0][1] - intPntProj[1]);
								d = (linepProj[0][0] - linepProj[1][0]) * (linepProj[0][0] - linepProj[1][0]) +
									(linepProj[0][1] - linepProj[1][1]) * (linepProj[0][1] - linepProj[1][1]);

								d1 = sqrt(d1);
								d = sqrt(d);
								w2 = d1 / d;
								w1 = 1.0 - w2;
								intPnt[0] = w1 * linep[0][0] + w2 * linep[1][0];
								intPnt[1] = w1 * linep[0][1] + w2 * linep[1][1];
								intPnt[2] = w1 * linep[0][2] + w2 * linep[1][2];
							}
						}
					}
				}
				else
				{/* 面片退化了，不处理 */
				//assert(false);
				//printf("Error in lin_tri_intersect3d(...). Degenerate triangle.\n");
					isInt = 1;
					*intTyp = LTI_INTERSECT_DEG_FACE;
				}
			}
			else if ((i1Ort == 0.0 && i2Ort != 0.0) || (i1Ort != 0.0 && i2Ort == 0.0))
			{
				/* -------------------------------------------------------------------
				 * 直线只有一个端点在面上时，这个点为交点
				 * ----------------------------------------------------------------*/
				intEndIdx = i1Ort == 0.0 ? 0 : 1;
				//不对，intTyp返回3，但其实此时相交的类型是点
				isInt = pnt_tri_intersect3d(linep[intEndIdx], facep, intTyp, intCod, intPnt);
			}
			else
			{
				/* -------------------------------------------------------------------
				 * 直线2个点都不在平面时，假设(p1,p2,p3)法向指向p4(直线的1个端点，另一个
				 * 负法向端点记为p5，则p4/p5和p1 p2 p3 相交的充要条件是：
				 * p1p2p5p4/p2p3p5p4/p3p1p5p4的体积都大于等于0 (相应orient3d则小于0)
				 * ----------------------------------------------------------------*/
				p4Idx = i1Ort > 0.0 ? 1 : 0;
				p5Idx = i1Ort > 0.0 ? 0 : 1;

				j1Ort = orient3d(facep[0], facep[1], linep[p4Idx], linep[p5Idx]);
				j2Ort = orient3d(facep[1], facep[2], linep[p4Idx], linep[p5Idx]);
				j3Ort = orient3d(facep[2], facep[0], linep[p4Idx], linep[p5Idx]);

				if (bEpsilon)
				{
					if (fabs(j1Ort) < 1.0e-15)
						j1Ort = 0.0;
					if (fabs(j2Ort) < 1.0e-15)
						j2Ort = 0.0;
					if (fabs(j3Ort) < 1.0e-15)
						j3Ort = 0.0;
				}

				/* 如果j1Ort/j2Ort/j3Ort有2个等于0 (因为2个点都不在面上，且面不退化，可排除3个都等于0的可能)，则交点在点上 */
				if (j1Ort == 0.0 && j2Ort == 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 1;
				}
				else if (j2Ort == 0.0 && j3Ort == 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 2;
				}
				else if (j3Ort == 0.0 && j1Ort == 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 0;
				}

				/* 如果相交，保存交点 */
				if (isInt)
				{
					if (intPnt)
					{
						intPnt[0] = facep[*intCod][0];
						intPnt[1] = facep[*intCod][1];
						intPnt[2] = facep[*intCod][2];
					}

					return isInt;
				}

				/* 如果j1Ort/j2Ort/j3Ort只有一个等于0，且另外两个都大于0，则交点在边上 */
				if (j1Ort == 0.0 && j2Ort > 0.0 && j3Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_EDG;
					*intCod = 0;
				}
				else if (j2Ort == 0.0 && j3Ort > 0.0 && j1Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_EDG;
					*intCod = 1;
				}
				else if (j3Ort == 0.0 && j1Ort > 0.0 && j2Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_EDG;
					*intCod = 2;
				}
				else if (j1Ort > 0.0 && j2Ort > 0.0 && j3Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_FAC;
				}

				/* 如果相交，求出交点 */
				if (isInt)
				{
					/* 当我们需要更精确的坐标时，我们需要更小心处理浮点误差 */
#if 0
					w1 = fabs(i2Ort) / (fabs(i1Ort) + fabs(i2Ort));
					w2 = 1.0 - w1;

					intPnt[0] = w1 * linep[0][0] + w2 * linep[1][0];
					intPnt[1] = w1 * linep[0][1] + w2 * linep[1][1];
					intPnt[2] = w1 * linep[0][2] + w2 * linep[1][2];
#else
					if (intPnt)
					{
						i1Ort = orient3d(facep[0], facep[1], facep[2], linep[0]);
						i2Ort = orient3d(facep[0], facep[1], facep[2], linep[1]);
						if (i1Ort < 0.0)
							i1Ort = -i1Ort;
						else
							i2Ort = -i2Ort;

						intPnt[0] = fixedSplitPoint(i1Ort, i2Ort, linep[0][0], linep[1][0]);
						intPnt[1] = fixedSplitPoint(i1Ort, i2Ort, linep[0][1], linep[1][1]);
						intPnt[2] = fixedSplitPoint(i1Ort, i2Ort, linep[0][2], linep[1][2]);
					}
#endif
				}
			}

			return isInt;
		}

		int lin_tri_intersect3d_idx(int iline[2], int iface[3], double linep[2][3], double facep[3][3],
			int* intTyp, int* intCod, double intPnt[3], bool bEpsilon)
		{
			int isInt = 0;
			bool ln1Similar, ln2Similar;
			double linepj[3][2] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };	/* 共面情形的投影点 */
			int m;
			int fcSimilarCnt = 0, iSimilarIdx;
			int iTemp;
			double dTemp;

			/* 面上有几个节点是共面的 */
			for (m = 0, fcSimilarCnt = 0; m < 3; m++)
			{
				if (iface[m] == iline[0] || iface[m] == iline[1])
				{
					iSimilarIdx = m;
					fcSimilarCnt++;
				}
			}

			switch (fcSimilarCnt)
			{
			case 0:
				isInt = lin_tri_intersect3d(linep, facep, intTyp, intCod, intPnt, bEpsilon);
				break;
			case 1:
				ln1Similar = (iline[0] == iface[0] ||
					iline[0] == iface[1] ||
					iline[0] == iface[2]);
				ln2Similar = (iline[1] == iface[0] ||
					iline[1] == iface[1] ||
					iline[1] == iface[2]);

				assert(ln1Similar || ln2Similar);
				if (ln2Similar)
				{
					assert(!ln1Similar);
					iTemp = iline[0];
					iline[0] = iline[1];
					iline[1] = iTemp;
					for (m = 0; m < 3; m++)
					{
						dTemp = linep[0][m];
						linep[0][m] = linep[1][m];
						linep[1][m] = dTemp;
					}
				}
				else
					assert(!ln2Similar);

				isInt = isintersect_oneSharePoint(iface, iline, facep, linep, intTyp, intCod, intPnt);

				break;
			case 2:
				*intTyp = LTI_INTERSECT_NUL;
				isInt = 0; /* 边为面片的子线段, 我们认为不相交 */
				break;
			case 3:
				//		printf("Error in lin_tri_intersect3d_idx(...). Degenerate triangle.\n");
				isInt = 1; /* 3点共线的面片，自然不允许 */
				*intTyp = LTI_INTERSECT_DEG_FACE;
				break;
			}

			return isInt;
		}

		//检查边-边相交，返回交点类型intTyp分别为交点在两个边上的情况，如果交点在端点上，intIdx返回交点点在iline中的索引，by zengyinjia
		//该版本为noOverlapPnt，即不检查两点索引不同但实际坐标相同的情况，可以加快速度，也不求出交点的具体位置，也可以加快速度。
		int lin_lin_intersect3d_idx_noOverlapPnt(int iline[2][2], double linep[2][2][3], int intTyp[2], int intIdx[2]) {
			//检查边退化成点的情况
			for (int i = 0; i < 2; i++) {
				if (iline[i][0] == iline[i][1]) {
					intTyp[0] = LTI_INTERSECT_NUL;
					intTyp[1] = LTI_INTERSECT_NUL;
					return 0;
				}
			}
			//检查是否共点
			constexpr int NodePair[4][2]{ {0,0},{0,1},{1,0},{1,1} };
			int nSamePoint = 0;
			for (int i = 0; i < 4; i++) {
				if (iline[0][NodePair[i][0]] == iline[1][NodePair[i][1]]) {
					nSamePoint++;
					//将相同点放到第一位
					std::swap(iline[0][NodePair[i][0]], iline[0][0]);
					std::swap(iline[1][NodePair[i][1]], iline[1][0]);
					for (int dim = 0; dim < 3; dim++) {
						//交换点坐标
						std::swap(linep[0][NodePair[i][0]][dim], linep[0][0][dim]);
						std::swap(linep[1][NodePair[i][1]][dim], linep[1][0][dim]);
					}
					break;
				}
			}
			//至此，若有相同的点，一定放在第一位上，因此如果第二位点相同，即为相同边（共享端点不算交点）
			if (iline[0][1] == iline[1][1]) {
				intTyp[0] = LTI_INTERSECT_NUL;
				intTyp[1] = LTI_INTERSECT_NUL;
				return 0;
			}

			//三个辅助点与linep[0][0][0]构成正交的向量基, TODO0 用非平凡的辅助向量进行优化
			double helperPnts[3][3] = { { linep[0][0][0] + 1, linep[0][0][1], linep[0][0][2] },
										{ linep[0][0][0], linep[0][0][1] + 1, linep[0][0][2] },
										{ linep[0][0][0], linep[0][0][1], linep[0][0][2] + 1 } };
			int ihelper;

			auto isColliear = [&helperPnts](double p1[3], double p2[3], double p3[3], int* ihelper) {
				//这里简单用叉乘求面积可能不准确，但是实现exact的版本又麻烦，因此通过辅助点的方式借用orient3d
				for (int i = 0; i < 3; i++) {
					if (orient3d(p1, p2, p3, helperPnts[i]) != 0) {
						*ihelper = i;
						return false;
					}
				}
				return true;
				};

			//相交检测
			if (nSamePoint) {
				//共享一个端点，这种情况下，只需检测是否共线
				if (isColliear(linep[0][0], linep[0][1], linep[1][1], &ihelper))
				{
					//共线
					double lineVec[2][3];
					for (int i = 0; i < 2; i++) {
						for (int dim = 0; dim < 3; dim++) {
							lineVec[i][dim] = linep[i][1][dim] - linep[i][0][dim];
						}
					}
					//点乘，判断是否共向
					double dot = 0;
					for (int dim = 0; dim < 3; dim++) {
						dot += lineVec[0][dim] * lineVec[1][dim];
					}
					if (dot <= 0) {
						//两边反向或其中一边退化为点，该情况认为没有交点
						intTyp[0] = LTI_INTERSECT_NUL;
						intTyp[1] = LTI_INTERSECT_NUL;
						return 0;
					}
					//判断长度，短边的端点1则为交点
					double len2[2]{ 0,0 };
					for (int dim = 0; dim < 3; dim++) {
						len2[0] += lineVec[0][dim] * lineVec[0][dim];
						len2[1] += lineVec[1][dim] * lineVec[1][dim];
					}
					//由于本函数假设不同索引点的坐标不同，所以两条边的第二个端点必不可能重合
					int ishorter = len2[0] > len2[1] ? 1 : 0;
					int ilonger = 1 - ishorter;
					intTyp[ilonger] = LTI_INTERSECT_EDG; //对长边来说，是在边中间相交的
					intTyp[ishorter] = LTI_INTERSECT_NOD; //对于短边来说，是在另外一个端点上相交的
					intIdx[ishorter] = iline[ishorter][1];
					return 1;
				}
				else {
					//不共线，肯定没有交点（共享端点不算交点）
					intTyp[0] = LTI_INTERSECT_NUL;
					intTyp[1] = LTI_INTERSECT_NUL;
					return 0;
				}
			}
			else {
				//没有共享端点
				//先判断是否共面
				double ori3d = orient3d(linep[0][0], linep[0][1], linep[1][0], linep[1][1]);
				if (ori3d != 0.0) {
					//若不共面，必不可能相交
					intTyp[0] = LTI_INTERSECT_NUL;
					intTyp[1] = LTI_INTERSECT_NUL;
					return 0;
				}

				//先判断是否共线
				if (isColliear(linep[0][0], linep[0][1], linep[1][0], &ihelper)) {
					if (isColliear(linep[0][0], linep[0][1], linep[1][1], &ihelper)) {
						//四点共线，检查一字型相交
						double lineVec[2][3];
						double lineLength2[2]{ 0,0 };
						for (int i = 0; i < 2; i++) {
							for (int dim = 0; dim < 3; dim++) {
								lineVec[i][dim] = linep[0][1][dim] - linep[0][0][dim];
								lineLength2[i] += lineVec[i][dim] * lineVec[i][dim];
							}
						}

						auto isPntOnLine = [](double linep1[3], double lineVec[3], double lineLength2, double pnt[3]) {
							double dot = 0;
							for (int dim = 0; dim < 3; dim++) {
								dot += lineVec[dim] * (pnt[dim] - linep1[dim]);
							}
							return dot > 0 && dot < lineLength2;
							};

						for (int i = 0; i < 2; i++) {
							for (int j = 0; j < 2; j++) {
								//判断边(1-i)的第j个点是否落在边i上
								if (isPntOnLine(linep[i][0], lineVec[i], lineLength2[i], linep[1 - i][j])) {
									intTyp[i] = LTI_INTERSECT_EDG;
									intTyp[1 - i] = LTI_INTERSECT_NOD;
									intIdx[1 - i] = iline[1 - i][j];
									//两边重叠的情况下，交点数量不好定义，这里就返回其中一个交点
									return 1;
								}
							}
						}
						//共线但不相交
						return 0;
					}
				}
				//若到此处没有返回，说明不四点共线，且ihelper存放的辅助点在四点所在平面之外

				double ori3dWithHelper[2][2]; //ori3dWithHelper[i][j]->边i两端点与边(1-i)的第j个端点，关于辅助点的orient3d值
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						ori3dWithHelper[i][j] = orient3d(linep[i][0], linep[i][1], linep[1 - i][j], helperPnts[ihelper]);
					}
				}

				// 检查T型相交
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						//此处i代表T一横的边索引，j代表一竖的边的端点索引
						//相交条件：点j在边i所在直线上，且边i的两端点在边(1-i)的两侧
						if (ori3dWithHelper[i][j] == 0.0 && ori3dWithHelper[1 - i][0] * ori3dWithHelper[1 - i][1] < 0.0) {
							intTyp[i] = LTI_INTERSECT_EDG;
							intTyp[1 - i] = LTI_INTERSECT_NOD;
							intIdx[1 - i] = iline[1 - i][j];
							return 1;
						}
					}
				}

				// 检查X型相交
				if (ori3dWithHelper[0][0] * ori3dWithHelper[0][1] < 0 && ori3dWithHelper[1][0] * ori3dWithHelper[1][1] < 0) {
					intTyp[0] = LTI_INTERSECT_EDG;
					intTyp[1] = LTI_INTERSECT_EDG;
					return 1;
				}

				// 不相交
				intTyp[0] = LTI_INTERSECT_NUL;
				intTyp[1] = LTI_INTERSECT_NUL;
				return 0;
			}
		}

		int lin_lin_intersect3d_exact(
			double linep[2][2][3],
			int intTyp[2],
			int intIdx[2],
			double intPnt[3]
		) {
			double* A = linep[0][0];
			double* B = linep[0][1];
			double* C = linep[1][0];
			double* D = linep[1][1];

			intTyp[0] = LTI_INTERSECT_NUL;
			intTyp[1] = LTI_INTERSECT_NUL;
			intIdx[0] = -1;
			intIdx[1] = -1;

			auto samePoint = [](double* P, double* Q) -> bool {
				return P[0] == Q[0] && P[1] == Q[1] && P[2] == Q[2];
				};

			auto copyPoint = [](double* dst, double* src) {
				dst[0] = src[0];
				dst[1] = src[1];
				dst[2] = src[2];
				};

			auto sqLen = [](double* P, double* Q) -> double {
				double x = P[0] - Q[0];
				double y = P[1] - Q[1];
				double z = P[2] - Q[2];
				return x * x + y * y + z * z;
				};

			// Degenerate segments.
			if (sqLen(A, B) == 0.0 || sqLen(C, D) == 0.0) {
				return 0;
			}

			// Endpoint sharing: unique intersection.
			if (samePoint(A, C)) {
				copyPoint(intPnt, A);
				intTyp[0] = LTI_INTERSECT_NOD;
				intTyp[1] = LTI_INTERSECT_NOD;
				intIdx[0] = 0;
				intIdx[1] = 0;
				return 1;
			}
			if (samePoint(A, D)) {
				copyPoint(intPnt, A);
				intTyp[0] = LTI_INTERSECT_NOD;
				intTyp[1] = LTI_INTERSECT_NOD;
				intIdx[0] = 0;
				intIdx[1] = 1;
				return 1;
			}
			if (samePoint(B, C)) {
				copyPoint(intPnt, B);
				intTyp[0] = LTI_INTERSECT_NOD;
				intTyp[1] = LTI_INTERSECT_NOD;
				intIdx[0] = 1;
				intIdx[1] = 0;
				return 1;
			}
			if (samePoint(B, D)) {
				copyPoint(intPnt, B);
				intTyp[0] = LTI_INTERSECT_NOD;
				intTyp[1] = LTI_INTERSECT_NOD;
				intIdx[0] = 1;
				intIdx[1] = 1;
				return 1;
			}

			// Check coplanarity exactly.
			// If four points are not coplanar, two 3D segments cannot intersect.
			double vol = orient3dexact(A, B, C, D);
			if (vol != 0.0) {
				return 0;
			}

			// Exact-style collinearity check using orient3dexact with auxiliary points.
			// Three collinear points will be coplanar with any auxiliary point.
			auto isCollinear3D = [](double* P, double* Q, double* R) -> bool {
				double H[3][3] = {
					{ P[0] + 1.0, P[1],       P[2]       },
					{ P[0],       P[1] + 1.0, P[2]       },
					{ P[0],       P[1],       P[2] + 1.0 }
				};

				for (int i = 0; i < 3; ++i) {
					if (orient3dexact(P, Q, R, H[i]) != 0.0) {
						return false;
					}
				}
				return true;
				};

			bool colC = isCollinear3D(A, B, C);
			bool colD = isCollinear3D(A, B, D);

			// Collinear case: intersection may be a point or an overlapping segment.
			if (colC && colD) {
				double dir[3] = {
					B[0] - A[0],
					B[1] - A[1],
					B[2] - A[2]
				};

				int dim = 0;
				if (fabs(dir[1]) > fabs(dir[dim])) dim = 1;
				if (fabs(dir[2]) > fabs(dir[dim])) dim = 2;

				double a0 = A[dim], a1 = B[dim];
				double c0 = C[dim], c1 = D[dim];

				if (a0 > a1) std::swap(a0, a1);
				if (c0 > c1) std::swap(c0, c1);

				double lo = std::max(a0, c0);
				double hi = std::min(a1, c1);

				if (lo > hi) {
					return 0;
				}

				// Touch at one point.
				if (lo == hi) {
					double* P = NULL;

					if (A[dim] == lo) P = A;
					else if (B[dim] == lo) P = B;
					else if (C[dim] == lo) P = C;
					else if (D[dim] == lo) P = D;

					if (P == NULL) {
						return 0;
					}

					copyPoint(intPnt, P);

					if (samePoint(P, A)) {
						intTyp[0] = LTI_INTERSECT_NOD;
						intIdx[0] = 0;
					}
					else if (samePoint(P, B)) {
						intTyp[0] = LTI_INTERSECT_NOD;
						intIdx[0] = 1;
					}
					else {
						intTyp[0] = LTI_INTERSECT_EDG;
						intIdx[0] = -1;
					}

					if (samePoint(P, C)) {
						intTyp[1] = LTI_INTERSECT_NOD;
						intIdx[1] = 0;
					}
					else if (samePoint(P, D)) {
						intTyp[1] = LTI_INTERSECT_NOD;
						intIdx[1] = 1;
					}
					else {
						intTyp[1] = LTI_INTERSECT_EDG;
						intIdx[1] = -1;
					}

					return 1;
				}

				// Overlap in a segment, not a unique point.
				intTyp[0] = LTI_INTERSECT_INS;
				intTyp[1] = LTI_INTERSECT_INS;
				return 2;
			}

			// Coplanar but not collinear.
			// Choose a 2D projection where AB and C are non-collinear.
			int id0 = -1, id1 = -1;
			const int proj[3][2] = {
				{0, 1},
				{1, 2},
				{2, 0}
			};

			double A2[2], B2[2], C2[2], D2[2];

			for (int k = 0; k < 3; ++k) {
				int i = proj[k][0];
				int j = proj[k][1];

				A2[0] = A[i]; A2[1] = A[j];
				B2[0] = B[i]; B2[1] = B[j];
				C2[0] = C[i]; C2[1] = C[j];

				if (orient2d(A2, B2, C2) != 0.0) {
					id0 = i;
					id1 = j;
					break;
				}
			}

			if (id0 == -1) {
				// Should not happen if not collinear, but keep safe.
				return 0;
			}

			A2[0] = A[id0]; A2[1] = A[id1];
			B2[0] = B[id0]; B2[1] = B[id1];
			C2[0] = C[id0]; C2[1] = C[id1];
			D2[0] = D[id0]; D2[1] = D[id1];

			double o1 = orient2d(A2, B2, C2);
			double o2 = orient2d(A2, B2, D2);
			double o3 = orient2d(C2, D2, A2);
			double o4 = orient2d(C2, D2, B2);

			// Segment intersection test in projected plane.
			if (o1 * o2 > 0.0 || o3 * o4 > 0.0) {
				return 0;
			}

			// If intersection is A or B.
			if (o3 == 0.0) {
				copyPoint(intPnt, A);
				intTyp[0] = LTI_INTERSECT_NOD;
				intIdx[0] = 0;

				if (samePoint(A, C)) {
					intTyp[1] = LTI_INTERSECT_NOD;
					intIdx[1] = 0;
				}
				else if (samePoint(A, D)) {
					intTyp[1] = LTI_INTERSECT_NOD;
					intIdx[1] = 1;
				}
				else {
					intTyp[1] = LTI_INTERSECT_EDG;
					intIdx[1] = -1;
				}

				return 1;
			}

			if (o4 == 0.0) {
				copyPoint(intPnt, B);
				intTyp[0] = LTI_INTERSECT_NOD;
				intIdx[0] = 1;

				if (samePoint(B, C)) {
					intTyp[1] = LTI_INTERSECT_NOD;
					intIdx[1] = 0;
				}
				else if (samePoint(B, D)) {
					intTyp[1] = LTI_INTERSECT_NOD;
					intIdx[1] = 1;
				}
				else {
					intTyp[1] = LTI_INTERSECT_EDG;
					intIdx[1] = -1;
				}

				return 1;
			}

			// Proper intersection or C/D lies on AB.
			// Use exact oriented areas to interpolate along AB.
			double sA = o3;
			double sB = o4;

			if (sA < 0.0) {
				sA = -sA;
			}
			else {
				sB = -sB;
			}

			intPnt[0] = fixedSplitPoint(sA, sB, A[0], B[0]);
			intPnt[1] = fixedSplitPoint(sA, sB, A[1], B[1]);
			intPnt[2] = fixedSplitPoint(sA, sB, A[2], B[2]);

			// Classify position on segment AB.
			if (samePoint(intPnt, A)) {
				intTyp[0] = LTI_INTERSECT_NOD;
				intIdx[0] = 0;
			}
			else if (samePoint(intPnt, B)) {
				intTyp[0] = LTI_INTERSECT_NOD;
				intIdx[0] = 1;
			}
			else {
				intTyp[0] = LTI_INTERSECT_EDG;
				intIdx[0] = -1;
			}

			// Classify position on segment CD.
			if (samePoint(intPnt, C)) {
				intTyp[1] = LTI_INTERSECT_NOD;
				intIdx[1] = 0;
			}
			else if (samePoint(intPnt, D)) {
				intTyp[1] = LTI_INTERSECT_NOD;
				intIdx[1] = 1;
			}
			else {
				intTyp[1] = LTI_INTERSECT_EDG;
				intIdx[1] = -1;
			}

			return 1;
		}

		inline int isintersect_oneSharePoint(int iface[], int iline[], double facep[][3], double linep[][3], int* intTyp, int* intCod, double intPnt[3])
		{
			int isInt = 0;
			double s1;
			int m;
			double facepProj[3][2], linepProj[2][2], intPntProj[2] = { 0 };
			double d, d1, w1, w2;

			/* 不共面，必不相交 */
			s1 = GEOM_FUNC::orient3d(facep[0], facep[1], facep[2], linep[1]);
			if (s1 == 0.0)
			{/* 共面? */
				/* -------------------------------------------------------------------
				 * 直线和平面共面时，情况很复杂，直线和面的边界可能有2个交点，此时
				 * 暂时只返回离直线首点近的交点
				 * ----------------------------------------------------------------*/
				if (proj_coplanr_points(facep[0], facep[1], facep[2], linep,
					facepProj[0], facepProj[1], facepProj[2], linepProj, 2) == 1)
				{
					isInt = lin_tri_intersect2d(linepProj, facepProj, intTyp, intCod, intPnt ? intPntProj : NULL);

					/* 将交点投影到平面上 */
					if (isInt == 1)
					{
						if (*intTyp == LTI_INTERSECT_NOD && intPnt)
						{/* 需要求交点位置 */
							intPnt[0] = facep[*intCod][0];
							intPnt[1] = facep[*intCod][1];
							intPnt[2] = facep[*intCod][2];
						}
						else if (isInt == 1 && intPnt)
						{
							d1 = (linepProj[0][0] - intPntProj[0]) * (linepProj[0][0] - intPntProj[0]) +
								(linepProj[0][1] - intPntProj[1]) * (linepProj[0][1] - intPntProj[1]);
							d = (linepProj[0][0] - linepProj[1][0]) * (linepProj[0][0] - linepProj[1][0]) +
								(linepProj[0][1] - linepProj[1][1]) * (linepProj[0][1] - linepProj[1][1]);

							d1 = sqrt(d1);
							d = sqrt(d);
							w2 = d1 / d;
							w1 = 1.0 - w2;
							intPnt[0] = w1 * linep[0][0] + w2 * linep[1][0];
							intPnt[1] = w1 * linep[0][1] + w2 * linep[1][1];
							intPnt[2] = w1 * linep[0][2] + w2 * linep[1][2];
						}
					}
				}
				else
				{/* 面片退化了，不处理 */
			//		printf("Error in isintersect_oneSharePoint(...). Degenerate triangle.\n");
					*intTyp = LTI_INTERSECT_DEG_FACE;
					isInt = 1;
				}
			}
			else
			{
				for (m = 0; m < 3; m++)
					if (iface[m] == iline[0] || iface[m] == iline[1])
						break;
				assert(m < 3);

				isInt = 1;
				*intTyp = NOD;
				*intCod = m;
				if (intPnt)
				{
					intPnt[0] = facep[m][0];
					intPnt[1] = facep[m][1];
					intPnt[2] = facep[m][2];
				}
			}

			return isInt;
		}

		/* -----------------------------------------------------------------------------------
		 * 这是一个判定2个三角面片是否相交的代码，它将在small polyhedron reconnection算法中调用
		 * -----------------------------------------------------------------------------------*/
		int tri_tri_intersect3d(int facei1[3], int facei2[3], double facep1[3][3], double facep2[3][3])
		{
			int sameVtxCnt = 0, sameVtx1[3], sameVtx2[3];
			int i, j;
			int isInt = 0;
			int ii1, ii2, ii3, ii4, ii5;
			double p1[3], p2[3], p3[3], p4[3];
			double dq1, dr1, dq2, dr2;

			/* -----------------------------------------------------------------
			 * 注意:
			 * 当pa/pb/pc的右手绕向指向pd时，调用orient3d(pa,pa,pc,pd)返回负数
			 *-----------------------------------------------------------------*/
			for (i = 0, sameVtxCnt = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					if (facei1[i] == facei2[j])
					{
						sameVtx1[sameVtxCnt] = i;
						sameVtx2[sameVtxCnt] = j;
						sameVtxCnt++;
						break;
					}
				}
			}

			switch (sameVtxCnt)
			{
			case 0:
				isInt = tri_tri_overlap_test_3d(facep1[0], facep1[1], facep1[2],
					facep2[0], facep2[1], facep2[2]);
				break;
			case 1:
				/* 假设共同的顶点为p，不同的顶点为q1/q2和r1/r2，则如果另外2个点在平面的同一侧，则不相交 */
				ii1 = sameVtx1[0];
				ii2 = (ii1 + 1) % 3;
				ii3 = (ii1 + 2) % 3;
				dq1 = orient3d(facep2[0], facep2[1], facep2[2], facep1[ii2]);
				dr1 = orient3d(facep2[0], facep2[1], facep2[2], facep1[ii3]);
				if (dq1 * dr1 > 0.0)
					break;	/* isInt = 0, 已经得到结果 */

				ii1 = sameVtx2[0];
				ii4 = (ii1 + 1) % 3;
				ii5 = (ii1 + 2) % 3;
				dq2 = orient3d(facep1[0], facep1[1], facep1[2], facep2[ii4]);
				dr2 = orient3d(facep1[0], facep1[1], facep1[2], facep2[ii5]);
				if (dq2 * dr2 > 0.0)
					break;	/* isInt = 0, 已经得到结果 */

				/* 如果q1或r1在face2的平面上，只需确认这个点和p1的联系不在其内角q2pr2内部 */
				if (dq1 * dr1 == 0.0)
				{
					/* 确定p1->p4不在内角p2p1p3内部 */
					isInt = line_cross_tri_coplane3d(facep2[ii1], facep2[ii4], facep2[ii5],
						dq1 == 0.0 ? facep1[ii2] : facep1[ii3]);
					if (dq1 == 0.0 && dr1 == 0.0 && isInt != 1)
					{
						isInt = line_cross_tri_coplane3d(facep2[ii1], facep2[ii4], facep2[ii5], facep1[ii3]);
						if (isInt == 1)	/* 等于0则还要继续 */
							break;
					}
					else
						break; /* 不共面时，已得出结论 */
				}
				assert(isInt == 0); /* 此时isInt必定等于0 */

				/* 同理，如果q2或r2在face1的平面上，只需确认这个点和p2的连线不在其内角q1pr1内部 */
				if (dq2 * dr2 == 0.0)// && dq2 + dr2 != 0.0)
				{
					/* 确定p1->p4不在内角p2p1p3内部 */
					isInt = line_cross_tri_coplane3d(facep1[sameVtx1[0]], facep1[ii2], facep1[ii3],
						dq2 == 0.0 ? facep2[ii4] : facep2[ii5]);
					if (dq2 == 0.0 && dr2 == 0.0 && isInt != 1)
					{
						isInt = line_cross_tri_coplane3d(facep1[sameVtx1[0]], facep1[ii2], facep1[ii3], facep2[ii5]);
					}
					break; /* 不管共面还是不共面，都已得到结论 */
				}

				/* 程序若进行到这，表明dp1*dq1 < 0.0 && dp2*dq2 < 0.0，这是只有判定
				 * q1(q2)->r1(r2)是否和f1(f2)相交
				 */
				assert(dq1 * dr1 < 0.0 && dq2 * dr2 < 0);
				isInt = orient3d(facep1[1], facep1[0], dq2 < 0.0 ? facep2[ii4] : facep2[ii5], dq2 < 0.0 ? facep2[ii5] : facep2[ii4]) < 0.0 &&
					orient3d(facep1[2], facep1[1], dq2 < 0.0 ? facep2[ii4] : facep2[ii5], dq2 < 0.0 ? facep2[ii5] : facep2[ii4]) < 0.0 &&
					orient3d(facep1[0], facep1[2], dq2 < 0.0 ? facep2[ii4] : facep2[ii5], dq2 < 0.0 ? facep2[ii5] : facep2[ii4]) < 0.0;
				if (isInt == 1)
					break; /* 如果相交，则退出 */

				isInt = orient3d(facep2[1], facep2[0], dq1 < 0.0 ? facep1[ii2] : facep1[ii3], dq1 < 0.0 ? facep1[ii3] : facep1[ii2]) < 0.0 &&
					orient3d(facep2[2], facep2[1], dq1 < 0.0 ? facep1[ii2] : facep1[ii3], dq1 < 0.0 ? facep1[ii3] : facep1[ii2]) < 0.0 &&
					orient3d(facep2[0], facep2[2], dq1 < 0.0 ? facep1[ii2] : facep1[ii3], dq1 < 0.0 ? facep1[ii3] : facep1[ii2]) < 0.0;
				break;
			case 2:
				/* 实质只有4个不同的点，如果这4个点不共面，则必不相交 */
				ii4 = 3 - sameVtx2[0] - sameVtx2[1];
				if (orient3d(facep1[0], facep1[1], facep1[2], facep2[ii4]) == 0.0)
				{
					/* 如果4个点共面，首先将它们向使得投影面积最大的方向投影，
					 * 设共享的两个点为p1/p2，另外2个点为p3和p4，则必须保证p3
					 * 和p4位于直线的两侧
					 * ---------------------------------------------------*/
					 /* 求得平面的法向，沿使得平面投影面积最大的方向投影平面和直线 */
					ii1 = 3 - sameVtx1[0] - sameVtx1[1];
					if (proj_four_coplanr_points(facep1[(ii1 + 1) % 3], facep1[(ii1 + 2) % 3], facep1[ii1], facep2[ii4], p1, p2, p3, p4) == 1)
					{
						/* 判断p4和p3是否位于p1和p2的两侧 */
						isInt = orient2d(p1, p2, p3) * orient2d(p1, p2, p4) >= 0.0;
					}
					else
					{/* 平面退化成一条直线 */
						isInt = 1; /* 在SPR中，退化成直线的平面总是不允许的 */
					}
				}
				break;
			}

			return isInt;
		}

		/* -----------------------------------------------------------------------------------
		 * 这是一个判定2个三角面片是否相交的代码，它将在small polyhedron reconnection算法中调用
		 * -----------------------------------------------------------------------------------*/
		int tri_tri_intersect3d_fast(int facei1[3], int facei2[3], double* facep1[3], double* facep2[3])
		{
			int sameVtxCnt = 0, sameVtx1[3], sameVtx2[3];
			int i, j;
			int isInt = 0;
			int ii1, ii2, ii3, ii4, ii5;
			double p1[3], p2[3], p3[3], p4[3];
			double dq1, dr1, dq2, dr2;

			/* -----------------------------------------------------------------
			 * 注意:
			 * 当pa/pb/pc的右手绕向指向pd时，调用orient3d(pa,pa,pc,pd)返回负数
			 *-----------------------------------------------------------------*/
			for (i = 0, sameVtxCnt = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					if (facei1[i] == facei2[j])
					{
						sameVtx1[sameVtxCnt] = i;
						sameVtx2[sameVtxCnt] = j;
						sameVtxCnt++;
						break;
					}
				}
			}

			switch (sameVtxCnt)
			{
			case 0:
				isInt = tri_tri_overlap_test_3d(facep1[0], facep1[1], facep1[2],
					facep2[0], facep2[1], facep2[2]);
				break;
			case 1:
				/* 假设共同的顶点为p，不同的顶点为q1/q2和r1/r2，则如果另外2个点在平面的同一侧，则不相交 */
				ii1 = sameVtx1[0];
				ii2 = (ii1 + 1) % 3;
				ii3 = (ii1 + 2) % 3;
				dq1 = orient3d(facep2[0], facep2[1], facep2[2], facep1[ii2]);
				dr1 = orient3d(facep2[0], facep2[1], facep2[2], facep1[ii3]);
				if (dq1 * dr1 > 0.0)
					break;	/* isInt = 0, 已经得到结果 */

				ii1 = sameVtx2[0];
				ii4 = (ii1 + 1) % 3;
				ii5 = (ii1 + 2) % 3;
				dq2 = orient3d(facep1[0], facep1[1], facep1[2], facep2[ii4]);
				dr2 = orient3d(facep1[0], facep1[1], facep1[2], facep2[ii5]);
				if (dq2 * dr2 > 0.0)
					break;	/* isInt = 0, 已经得到结果 */

				/* 如果q1或r1在face2的平面上，只需确认这个点和p1的联系不在其内角q2pr2内部 */
				if (dq1 * dr1 == 0.0)
				{
					/* 确定p1->p4不在内角p2p1p3内部 */
					isInt = line_cross_tri_coplane3d(facep2[ii1], facep2[ii4], facep2[ii5],
						dq1 == 0.0 ? facep1[ii2] : facep1[ii3]);
					if (dq1 == 0.0 && dr1 == 0.0 && isInt != 1)
					{
						isInt = line_cross_tri_coplane3d(facep2[ii1], facep2[ii4], facep2[ii5], facep1[ii3]);
						if (isInt == 1)	/* 等于0则还要继续 */
							break;
					}
					else
						break; /* 不共面时，已得出结论 */
				}
				assert(isInt == 0); /* 此时isInt必定等于0 */

				/* 同理，如果q2或r2在face1的平面上，只需确认这个点和p2的连线不在其内角q1pr1内部 */
				if (dq2 * dr2 == 0.0)// && dq2 + dr2 != 0.0)
				{
					/* 确定p1->p4不在内角p2p1p3内部 */
					isInt = line_cross_tri_coplane3d(facep1[sameVtx1[0]], facep1[ii2], facep1[ii3],
						dq2 == 0.0 ? facep2[ii4] : facep2[ii5]);
					if (dq2 == 0.0 && dr2 == 0.0 && isInt != 1)
					{
						isInt = line_cross_tri_coplane3d(facep1[sameVtx1[0]], facep1[ii2], facep1[ii3], facep2[ii5]);
					}
					break; /* 不管共面还是不共面，都已得到结论 */
				}

				/* 程序若进行到这，表明dp1*dq1 < 0.0 && dp2*dq2 < 0.0，这是只有判定
				 * q1(q2)->r1(r2)是否和f1(f2)相交
				 */
				assert(dq1 * dr1 < 0.0 && dq2 * dr2 < 0);
				isInt = orient3d(facep1[1], facep1[0], dq2 < 0.0 ? facep2[ii4] : facep2[ii5], dq2 < 0.0 ? facep2[ii5] : facep2[ii4]) < 0.0 &&
					orient3d(facep1[2], facep1[1], dq2 < 0.0 ? facep2[ii4] : facep2[ii5], dq2 < 0.0 ? facep2[ii5] : facep2[ii4]) < 0.0 &&
					orient3d(facep1[0], facep1[2], dq2 < 0.0 ? facep2[ii4] : facep2[ii5], dq2 < 0.0 ? facep2[ii5] : facep2[ii4]) < 0.0;
				if (isInt == 1)
					break; /* 如果相交，则退出 */

				isInt = orient3d(facep2[1], facep2[0], dq1 < 0.0 ? facep1[ii2] : facep1[ii3], dq1 < 0.0 ? facep1[ii3] : facep1[ii2]) < 0.0 &&
					orient3d(facep2[2], facep2[1], dq1 < 0.0 ? facep1[ii2] : facep1[ii3], dq1 < 0.0 ? facep1[ii3] : facep1[ii2]) < 0.0 &&
					orient3d(facep2[0], facep2[2], dq1 < 0.0 ? facep1[ii2] : facep1[ii3], dq1 < 0.0 ? facep1[ii3] : facep1[ii2]) < 0.0;
				break;
			case 2:
				/* 实质只有4个不同的点，如果这4个点不共面，则必不相交 */
				ii4 = 3 - sameVtx2[0] - sameVtx2[1];
				if (orient3d(facep1[0], facep1[1], facep1[2], facep2[ii4]) == 0.0)
				{
					/* 如果4个点共面，首先将它们向使得投影面积最大的方向投影，
					 * 设共享的两个点为p1/p2，另外2个点为p3和p4，则必须保证p3
					 * 和p4位于直线的两侧
					 * ---------------------------------------------------*/
					 /* 求得平面的法向，沿使得平面投影面积最大的方向投影平面和直线 */
					ii1 = 3 - sameVtx1[0] - sameVtx1[1];
					if (proj_four_coplanr_points(facep1[(ii1 + 1) % 3], facep1[(ii1 + 2) % 3], facep1[ii1], facep2[ii4], p1, p2, p3, p4) == 1)
					{
						/* 判断p4和p3是否位于p1和p2的两侧 */
						isInt = orient2d(p1, p2, p3) * orient2d(p1, p2, p4) >= 0.0;
					}
					else
					{/* 平面退化成一条直线 */
						isInt = 1; /* 在SPR中，退化成直线的平面总是不允许的 */
					}
				}
				break;
			}

			return isInt;
		}

		double distance_square(double p1[3], double p2[3])
		{
			double dx = p2[0] - p1[0];
			double dy = p2[1] - p1[1];
			double dz = p2[2] - p1[2];
			return dx * dx + dy * dy + dz * dz;
		}

		/* 计算一个四面体单元的形状质量值 */
		extern double tetrahedron_gamma(double p1[3], double p2[3], double p3[3], double p4[3])
		{
			double v = orient3d(p4, p1, p2, p3);
			double AB, BC, AC, AD, BD, CD;
			const double SQRT_3_72 = 72.0 * sqrt(3.0);

			/* 计算质量，确定是否可行 */
			AB = distance_square(p1, p2);
			BC = distance_square(p2, p3);
			AC = distance_square(p1, p3);
			AD = distance_square(p1, p4);
			BD = distance_square(p2, p4);
			CD = distance_square(p3, p4);

			v /= pow(AB + BC + AC + AD + BD + CD, 1.5) * SQRT_3_72;//算出这个四面体的质量

			return v;
		}
	};

	namespace GEOM_FUNC {
		//#zyj 下面的三个函数是能返回所有的交点的特别版本，用于面网格修复模块
		int proj_coplanr_points(double p1[3], double p2[3], double p3[3], double p4[][3], double proj1[2], double proj2[2], double proj3[2], double proj4[][2], int nOth);
#define ORIENT3D_FUN orient3d
#define ORIENT2D_FUN orient2d
		/* -----------------------------------------------------------------------------------
			* 这是一个判定1条线段和三角形是否相交的代码(3D)，并返回交点位置（同时返回两个交点的版本）
			* 相交类型：
			* PNT 相交于1个点（intCod返回点的编号0~2，i代表面的第i个顶点）：
			* EDG 相交于1条边（intCod返回边的编号0~2，i代表(i,(i+1)%3形成的边)
			* FAC 相交于1个面
			* intPnt返回交点的值
			* -----------------------------------------------------------------------------------*/
		int lin_tri_intersect3d_V2(double linep[2][3], double facep[3][3], int* intCnt, int intTyp[2], int intCod[2], double intPnt[2][3], double tolerance)
		{
			auto oriIsZero = [tolerance](double orient) {
				return std::abs(orient) <= tolerance;
				};
			//#TEST
			auto oriIsZero_Coplane = [tolerance](double orient) {
				return std::abs(orient) <= tolerance * 1e2;
				};

			double i1Ort, i2Ort;
			int isInt = 0;
			int intEndIdx, p4Idx, p5Idx;
			double j1Ort, j2Ort, j3Ort;
			double d1, d, w1, w2;
			double linepProj[2][2], facepProj[3][2], intPntProj[2][2];
			*intCnt = -1;
			*intTyp = LTI_INTERSECT_NUL;
			*intCod = -1;

			i1Ort = ORIENT3D_FUN(facep[0], facep[1], facep[2], linep[0]);
			i2Ort = ORIENT3D_FUN(facep[0], facep[1], facep[2], linep[1]);

			/* ---------------------------------------------------------------------
			 * 直线的端点在面的两侧时，必不相交
			 * --------------------------------------------------------------------*/
			if (i1Ort * i2Ort > 0.0)
			{
				return 0;
			}

			//#TEST
			if (oriIsZero_Coplane(i1Ort) && oriIsZero_Coplane(i2Ort))
			{
				/* -------------------------------------------------------------------
				 * 直线和平面共面时，情况很复杂，直线和面的边界可能有2个交点，此时
				 * 暂时只返回离直线首点近的交点
				 * ----------------------------------------------------------------*/
				if (proj_coplanr_points(facep[0], facep[1], facep[2], linep,
					facepProj[0], facepProj[1], facepProj[2], linepProj, 2) == 1)
				{
					isInt = lin_tri_intersect2d_V2(linepProj, facepProj, intCnt, intTyp, intCod, intPntProj, tolerance);

					/* 将交点投影到平面上 */
					for (int iInt = 0; iInt < *intCnt; iInt++) {
						if (intTyp[iInt] == LTI_INTERSECT_NOD)
						{
							intPnt[iInt][0] = facep[intCod[iInt]][0];
							intPnt[iInt][1] = facep[intCod[iInt]][1];
							intPnt[iInt][2] = facep[intCod[iInt]][2];
						}
						else
						{
							/* 当我们需要更精确的坐标时，我们需要更小心处理浮点误差 */

							d1 = (linepProj[0][0] - intPntProj[iInt][0]) * (linepProj[0][0] - intPntProj[iInt][0]) +
								(linepProj[0][1] - intPntProj[iInt][1]) * (linepProj[0][1] - intPntProj[iInt][1]);
							d = (linepProj[0][0] - linepProj[1][0]) * (linepProj[0][0] - linepProj[1][0]) +
								(linepProj[0][1] - linepProj[1][1]) * (linepProj[0][1] - linepProj[1][1]);

							d1 = sqrt(d1);
							d = sqrt(d);
							w2 = d1 / d;
							w1 = 1.0 - w2;
							intPnt[iInt][0] = w1 * linep[0][0] + w2 * linep[1][0];
							intPnt[iInt][1] = w1 * linep[0][1] + w2 * linep[1][1];
							intPnt[iInt][2] = w1 * linep[0][2] + w2 * linep[1][2];
						}
					}
				}
				else
				{/* 面片退化了，不处理 */
		//			assert(false);
				//	printf("Error in lin_tri_intersect3d(...). Degenerate triangle.\n");
					isInt = 1;
					intCnt[0] = 1;
					intTyp[0] = LTI_INTERSECT_DEG_FACE;
				}
			}
			else if ((oriIsZero(i1Ort) && i2Ort != 0.0) || (i1Ort != 0.0 && oriIsZero(i2Ort)))
			{
				/* -------------------------------------------------------------------
				 * 直线只有一个端点在面上时，这个点为交点
				 * ----------------------------------------------------------------*/
				intEndIdx = oriIsZero(i1Ort) ? 0 : 1;
				isInt = pnt_tri_intersect3d_V2(linep[intEndIdx], facep, intTyp, intCod, intPnt[0], tolerance);
				*intCnt = isInt ? 1 : 0;
			}
			else
			{
				/* -------------------------------------------------------------------
				 * 直线2个点都不在平面时，假设(p1,p2,p3)法向指向p4(直线的1个端点，另一个
				 * 负法向端点记为p5，则p4/p5和p1 p2 p3 相交的充要条件是：
				 * p1p2p5p4/p2p3p5p4/p3p1p5p4的体积都大于等于0 (相应orient3d则小于0)
				 * ----------------------------------------------------------------*/
				p4Idx = i1Ort > 0.0 ? 1 : 0;
				p5Idx = i1Ort > 0.0 ? 0 : 1;

				j1Ort = ORIENT3D_FUN(facep[0], facep[1], linep[p4Idx], linep[p5Idx]);
				j2Ort = ORIENT3D_FUN(facep[1], facep[2], linep[p4Idx], linep[p5Idx]);
				j3Ort = ORIENT3D_FUN(facep[2], facep[0], linep[p4Idx], linep[p5Idx]);

				/* 如果j1Ort/j2Ort/j3Ort有2个等于0 (因为2个点都不在面上，且面不退化，可排除3个都等于0的可能)，则交点在点上 */
				if (oriIsZero(j1Ort) && oriIsZero(j2Ort))
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 1;
				}
				else if (oriIsZero(j2Ort) && oriIsZero(j3Ort))
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 2;
				}
				else if (oriIsZero(j3Ort) && oriIsZero(j1Ort))
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 0;
				}

				/* 如果相交，保存交点 */
				if (isInt)
				{
					*intCnt = 1;
					if (intPnt)
					{
						intPnt[0][0] = facep[*intCod][0];
						intPnt[0][1] = facep[*intCod][1];
						intPnt[0][2] = facep[*intCod][2];
					}

					return isInt;
				}

				/* 如果j1Ort/j2Ort/j3Ort只有一个等于0，且另外两个都大于0，则交点在边上 */
				if (oriIsZero(j1Ort) && j2Ort > 0.0 && j3Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_EDG;
					*intCod = 0;
				}
				else if (oriIsZero(j2Ort) && j3Ort > 0.0 && j1Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_EDG;
					*intCod = 1;
				}
				else if (oriIsZero(j3Ort) && j1Ort > 0.0 && j2Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_EDG;
					*intCod = 2;
				}
				else if (j1Ort > 0.0 && j2Ort > 0.0 && j3Ort > 0.0)
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_FAC;
				}

				/* 如果相交，求出交点 */
				if (isInt)
				{
					*intCnt = 1;
					/* 当我们需要更精确的坐标时，我们需要更小心处理浮点误差 */
#if 0
					w1 = fabs(i2Ort) / (fabs(i1Ort) + fabs(i2Ort));
					w2 = 1.0 - w1;

					intPnt[0] = w1 * linep[0][0] + w2 * linep[1][0];
					intPnt[1] = w1 * linep[0][1] + w2 * linep[1][1];
					intPnt[2] = w1 * linep[0][2] + w2 * linep[1][2];
#else
					if (intPnt)
					{
						i1Ort = ORIENT3D_FUN(facep[0], facep[1], facep[2], linep[0]);
						i2Ort = ORIENT3D_FUN(facep[0], facep[1], facep[2], linep[1]);
						if (i1Ort < 0.0)
							i1Ort = -i1Ort;
						else
							i2Ort = -i2Ort;

						intPnt[0][0] = fixedSplitPoint(i1Ort, i2Ort, linep[0][0], linep[1][0]);
						intPnt[0][1] = fixedSplitPoint(i1Ort, i2Ort, linep[0][1], linep[1][1]);
						intPnt[0][2] = fixedSplitPoint(i1Ort, i2Ort, linep[0][2], linep[1][2]);
					}
#endif
				}
			}

			return isInt;
		}

		/* -----------------------------------------------------------------------------------
		 * 这是一个判定1条线段和三角形是否相交的代码(2D)，并返回交点位置
		 * 相交类型：
		 * PNT 相交于1个点（intCod返回点的编号0~2，i代表面的第i个顶点）：
		 * EDG 相交于1条边（intCod返回边的编号0~2，i代表(i,(i+1)%3形成的边)
		 * FAC 相交于1个面
		 * intPnt返回交点的值
		 * 注意：在共面情形下，一条线可能会和一个面的两条边都相交，该版本支持两个交点的情况
		 * -----------------------------------------------------------------------------------*/
		int lin_tri_intersect2d_V2(double linep[2][2], double facep[3][2], int* intCnt, int intTyp[2], int intCod[2], double intPnt[2][2], double tolerance)
		{
			auto oriIsZero = [tolerance](double orient) {
				return std::abs(orient) <= tolerance;
				};
			auto modifyOri = [&oriIsZero](double orient) {
				if (oriIsZero(orient)) {
					return 0.0;
				}
				return orient;
				};

			int isInt = 0, isInner1, isInner2;
			double iOrt, iOrt1, iOrt2, iOrt3, iOrt11, iOrt12, iOrt13, iOrt21, iOrt22, iOrt23, kOrt1, kOrt2;
			double iOrtArr[3], jOrtArr[3], mOrtArr[3]; /* mOrtArr:面的3个顶点相对于直线的绕向 */
			double w1 = 0.0, w2 = 0.0;
			int m, nodCod, edgCod, candEdgNum = 0, instNodNum = 0;
			double candEdgCod[3], twoIntPnt[2][2], twoIntCod[2], twoIntTyp[2]; /* 可能有2个交点 */
			double d1 = 0.0, d2 = 0.0;
			int iCord;
			*intCnt = -1;
			/* -------------------------------------------------------------------
			 * 判定2个顶点是否在三角形内部
			 * ----------------------------------------------------------------*/
			iOrt = modifyOri(ORIENT2D_FUN(facep[0], facep[1], facep[2]));
			iOrt11 = modifyOri(ORIENT2D_FUN(facep[0], facep[1], linep[0]));
			iOrt12 = modifyOri(ORIENT2D_FUN(facep[1], facep[2], linep[0]));
			iOrt13 = modifyOri(ORIENT2D_FUN(facep[2], facep[0], linep[0]));
			iOrt21 = modifyOri(ORIENT2D_FUN(facep[0], facep[1], linep[1]));
			iOrt22 = modifyOri(ORIENT2D_FUN(facep[1], facep[2], linep[1]));
			iOrt23 = modifyOri(ORIENT2D_FUN(facep[2], facep[0], linep[1]));

			isInner1 = iOrt * iOrt11 >= 0 && iOrt * iOrt12 >= 0 && iOrt * iOrt13 >= 0;
			isInner2 = iOrt * iOrt21 >= 0 && iOrt * iOrt22 >= 0 && iOrt * iOrt23 >= 0;

			*intTyp = LTI_INTERSECT_NUL; /* 预设为不相交 */

			if (isInner1 && isInner2)
			{
				/* 直线两个点都在平面的内部或边界上时，这条直线被包含在平面内部 */
				isInt = 1;
				*intCnt = 2;

				double vOrt[2][3] = {
					{iOrt11, iOrt12, iOrt13},
					{iOrt21, iOrt22, iOrt23}
				};

				for (int iInt = 0; iInt < 2; iInt++) {
					int nZeroOri = 0;
					int iZeroOri, iNotZeroOri;
					for (int dim = 0; dim < 3; dim++) {
						if (oriIsZero(vOrt[iInt][dim])) {
							nZeroOri++;
							iZeroOri = dim;
						}
						else {
							iNotZeroOri = dim;
						}
					}
					if (nZeroOri == 0) { //在面内
						intTyp[iInt] = LTI_INTERSECT_FAC;
					}
					if (nZeroOri == 1) { //在边上
						intTyp[iInt] = LTI_INTERSECT_EDG;
						intCod[iInt] = iZeroOri;
					}
					if (nZeroOri == 2) { //在点上
						intTyp[iInt] = LTI_INTERSECT_NOD;
						intCod[iInt] = (iNotZeroOri + 2) % 3;
					}
					if (nZeroOri == 3) {
						throw std::logic_error("Self-intersections occur within tolerance. Fix failed!");
					}
					intPnt[iInt][0] = linep[iInt][0];
					intPnt[iInt][1] = linep[iInt][1];
				}
			}
			else if (isInner1 || isInner2)
			{
				/* 直线1个点在平面的内部或边界上时，这条直线和平面必有一个交点 */
				int innPt = isInner1 ? 0 : 1;
				iOrtArr[0] = iOrt1 = isInner1 ? iOrt11 : iOrt21;
				iOrtArr[1] = iOrt2 = isInner1 ? iOrt12 : iOrt22;
				iOrtArr[2] = iOrt3 = isInner1 ? iOrt13 : iOrt23;
				jOrtArr[0] = isInner2 ? iOrt11 : iOrt21;
				jOrtArr[1] = isInner2 ? iOrt12 : iOrt22;
				jOrtArr[2] = isInner2 ? iOrt13 : iOrt23;

				/* 如果iOrt1/iOrt2/iOrt3有2个等于0 (投影成功，则排除3个都等于0的可能)，则交点在点上 */
				if (oriIsZero(iOrt1) && oriIsZero(iOrt2))
				{
					assert(iOrt3 != 0.0);
					isInt = 1;
					intTyp[0] = LTI_INTERSECT_NOD;
					intCod[0] = 1;
					intPnt[0][0] = linep[innPt][0];
					intPnt[0][1] = linep[innPt][1];
				}
				else if (oriIsZero(iOrt2) && oriIsZero(iOrt3))
				{
					assert(iOrt1 != 0.0);
					isInt = 1;
					intTyp[0] = LTI_INTERSECT_NOD;
					intCod[0] = 2;
					intPnt[0][0] = linep[innPt][0];
					intPnt[0][1] = linep[innPt][1];
				}
				else if (oriIsZero(iOrt3) && oriIsZero(iOrt1))
				{
					assert(iOrt2 != 0.0);
					isInt = 1;
					intTyp[0] = LTI_INTERSECT_NOD;
					intCod[0] = 0;
					intPnt[0][0] = linep[innPt][0];
					intPnt[0][1] = linep[innPt][1];
				}

				/* 如果iOrt1/iOrt2/iOrt3只有一个等于0，且另外两个与iOrt同号，则交点在边上 */
				if (isInt == 0)
				{
					if (oriIsZero(iOrt1))
					{
						/* 走到这一步，iOrt2和iOrt3值必然和iOrt同号 */
						assert(iOrt2 * iOrt > 0.0 && iOrt3 * iOrt > 0.0);
						isInt = 1;
						intTyp[0] = LTI_INTERSECT_EDG;
						intCod[0] = 0;
						intPnt[0][0] = linep[innPt][0];
						intPnt[0][1] = linep[innPt][1];
					}
					if (oriIsZero(iOrt2))
					{
						/* 走到这一步，iOrt3和iOrt1值必然和iOrt同号 */
						assert(iOrt3 * iOrt > 0.0 && iOrt1 * iOrt > 0.0);
						isInt = 1;
						intTyp[0] = LTI_INTERSECT_EDG;
						intCod[0] = 1;
						intPnt[0][0] = linep[innPt][0];
						intPnt[0][1] = linep[innPt][1];
					}
					if (oriIsZero(iOrt3))
					{
						/* 走到这一步，iOrt1和iOrt2值必然和iOrt同号 */
						assert(iOrt1 * iOrt > 0.0 && iOrt2 * iOrt > 0.0);
						isInt = 1;
						intTyp[0] = LTI_INTERSECT_EDG;
						intCod[0] = 2;
						intPnt[0][0] = linep[innPt][0];
						intPnt[0][1] = linep[innPt][1];
					}
				}

				/* 内部点不在点上或边上，只有一个交点在边界 */
				if (isInt == 0)
				{
					*intCnt = 2;
					mOrtArr[0] = ORIENT2D_FUN(linep[0], linep[1], facep[0]);
					mOrtArr[1] = ORIENT2D_FUN(linep[0], linep[1], facep[1]);
					mOrtArr[2] = ORIENT2D_FUN(linep[0], linep[1], facep[2]);

					/* 第一类情形，面有1个端点在直线上 （走到这，面不可能有2个端点在直线上）*/
					for (nodCod = 0; nodCod < 3; nodCod++)
					{
						if (oriIsZero(mOrtArr[nodCod]))
						{
							isInt = 1;

							intTyp[0] = LTI_INTERSECT_NOD;
							intCod[0] = nodCod;
							intPnt[0][0] = facep[nodCod][0];
							intPnt[0][1] = facep[nodCod][1];

							break;
						}
					}

					/* 第二类情形，面有1条边和它相交 */
					if (isInt == 0)
					{
						for (edgCod = 0; edgCod < 3; edgCod++)
						{
							if (mOrtArr[edgCod] * mOrtArr[(edgCod + 1) % 3] < 0.0 && iOrtArr[edgCod] * jOrtArr[edgCod] < 0.0)
							{/* 相交 */
								isInt = 1;

								intTyp[0] = LTI_INTERSECT_EDG;
								intCod[0] = edgCod;

								w1 = isInner1 ? jOrtArr[edgCod] : iOrtArr[edgCod];
								w2 = isInner1 ? iOrtArr[edgCod] : jOrtArr[edgCod];
								w1 = fabs(w1) / (fabs(w1) + fabs(w2));
								w2 = 1.0 - w1;
								intPnt[0][0] = w1 * linep[0][0] + w2 * linep[1][0];
								intPnt[0][1] = w1 * linep[0][1] + w2 * linep[1][1];
								break;
							}
						}
					}

					intTyp[1] = LTI_INTERSECT_FAC;
					intPnt[1][0] = linep[innPt][0];
					intPnt[1][1] = linep[innPt][1];
				}

				/* 否则，可能有两个交点 */
				else
				{
					isInt = 1;
					/* 显然，直线可能还会和平面有第2个交点 */
					if (intTyp[0] == LTI_INTERSECT_NOD)
					{
						/* 判定直线是否在顶点所在的内角上 */
						assert(intCod[0] >= 0 && intCod[0] <= 2);
						nodCod = intCod[0];
						/* 在内角上 */
						if (jOrtArr[nodCod] * iOrt >= 0.0 && jOrtArr[(nodCod + 2) % 3] * iOrt >= 0.0)
						{
							/* 相交 */
							*intCnt = 2;
							if (oriIsZero(jOrtArr[nodCod]))
							{
								/* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
								/* 边的编号是edgCod, 两个顶点是edgCod和(edgCod+1)%3 */
								edgCod = nodCod;
								intTyp[1] = LTI_INTERSECT_NOD;
								intCod[1] = edgCod == nodCod ? (edgCod + 1) % 3 : edgCod;
								intPnt[1][0] = facep[intCod[1]][0];
								intPnt[1][1] = facep[intCod[1]][1];
							}
							else if (oriIsZero(jOrtArr[(nodCod + 2) % 3]))
							{/* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
								/* 边的编号是edgCod, 两个顶点是edgCod和(edgCod+1)%3 */
								edgCod = (nodCod + 2) % 3;
								intTyp[1] = LTI_INTERSECT_NOD;
								intCod[1] = edgCod == nodCod ? (edgCod + 1) % 3 : edgCod;
								intPnt[1][0] = facep[intCod[1]][0];
								intPnt[1][1] = facep[intCod[1]][1];
							}
							else
							{/* 交点在边上 */
								/* 当第1个点在平面内部，对应的jOrt代表不在平面点的和第3条边的面积，它越大，
								 * 则直线第1个点的贡献越大
								 */
								edgCod = (nodCod + 1) % 3;
								intTyp[1] = LTI_INTERSECT_EDG;
								intCod[1] = edgCod;

								w1 = isInner1 ? jOrtArr[edgCod] : iOrtArr[edgCod];
								w2 = isInner1 ? iOrtArr[edgCod] : jOrtArr[edgCod];
								w1 = fabs(w1) / (fabs(w1) + fabs(w2));
								w2 = 1.0 - w1;

								intPnt[1][0] = w1 * linep[0][0] + w2 * linep[1][0];
								intPnt[1][1] = w1 * linep[0][1] + w2 * linep[1][1];
							}
						}
						/* 不在内角上，一个交点 */
						else {
							*intCnt = 1;
						}
					}
					else if (*intTyp == LTI_INTERSECT_EDG)
					{
						/* 可能和两位两条边相交 */
						assert(*intCod >= 0 && *intCod <= 2);
						edgCod = intCod[0];
						/* 如果另外1个顶点在边的正侧，存在两个交点 */
						if (jOrtArr[edgCod] * iOrt > 0.0)
						{
							/* 如果另外2条边对应的jOrt值都小于0，还有进一步区分 */
							*intCnt = 2;
							if (jOrtArr[(edgCod + 1) % 3] * iOrt < 0.0 && jOrtArr[(edgCod + 2) % 3] * iOrt < 0.0)
							{
								/* 判定edgCod + 1的两个端点是否在直线的两侧 */
								edgCod = (edgCod + 1) % 3;
								kOrt1 = ORIENT2D_FUN(linep[0], linep[1], facep[edgCod]);
								kOrt2 = ORIENT2D_FUN(linep[0], linep[1], facep[(edgCod + 1) % 3]);
								if (oriIsZero(kOrt1) || oriIsZero(kOrt2))
								{
									/* 相交于面片的一个顶点 */
									isInt = 1;
									intTyp[1] = LTI_INTERSECT_NOD;
									intCod[1] = oriIsZero(kOrt1) ? edgCod : (edgCod + 1) % 3;
									intPnt[1][0] = facep[intCod[1]][0];
									intPnt[1][1] = facep[intCod[1]][1];
									edgCod = -1;
								}
								/* edgCod + 1的两个端点在直线的两侧，相交边为edgCod + 2 */
								else if (kOrt1 * kOrt2 > 0.0)
									edgCod = (edgCod + 2) % 3;
							}
							else if (jOrtArr[(edgCod + 1) % 3] * iOrt < 0.0)
							{
								edgCod = (edgCod + 1) % 3;
							}
							else
							{
								assert(jOrtArr[(edgCod + 2) % 3] * iOrt < 0.0);
								edgCod = (edgCod + 2) % 3;
							}

							if (edgCod >= 0)
							{
								/* 和平面第edgCod条边相交 */
								/* 当第1个点在平面内部，对应的jOrt代表不在平面点的和第3条边的面积，它越大，
								 * 则直线第1个点的贡献越大
								 */

								intTyp[1] = LTI_INTERSECT_EDG;
								intCod[1] = edgCod;

								w1 = isInner1 ? jOrtArr[edgCod] : iOrtArr[edgCod];
								w2 = isInner1 ? iOrtArr[edgCod] : jOrtArr[edgCod];
								w1 = fabs(w1) / (fabs(w1) + fabs(w2));
								w2 = 1.0 - w1;

								intPnt[1][0] = w1 * linep[0][0] + w2 * linep[1][0];
								intPnt[1][1] = w1 * linep[0][1] + w2 * linep[1][1];
							}
						}
						else if (oriIsZero(jOrtArr[edgCod]))
						{
							*intCnt = 2;
							/* 交点在顶点上，此时，三角形的一条边是原线段的子线段 */
							intTyp[1] = LTI_INTERSECT_NOD;
							/* 看下一条边的面积值，如果小于0，则是边edgCod的末点 */
							intCod[1] = jOrtArr[(edgCod + 1) % 3] < 0.0 ? (edgCod + 1) % 3 : edgCod;

							intPnt[1][0] = facep[intCod[1]][0];
							intPnt[1][1] = facep[intCod[1]][1];
						}
						/* 另外一个顶点在边的反侧，则只有一个顶点 */
						else {
							*intCnt = 1;
						}
					}
				}/* else (isInt != 0) */
			}/* else if (isInner1 || isInner2) */
			else {
				/* 直线两个点都在外部，仍然可能和平面相交于1个或2个点 */
				iOrtArr[0] = iOrt11;
				iOrtArr[1] = iOrt12;
				iOrtArr[2] = iOrt13;
				jOrtArr[0] = iOrt21;
				jOrtArr[1] = iOrt22;
				jOrtArr[2] = iOrt23;

				if (iOrtArr[0] * jOrtArr[0] <= 0.0 || iOrtArr[1] * jOrtArr[1] <= 0.0 || iOrtArr[2] * jOrtArr[2] <= 0.0)
				{/* 有可能相交 */
					/* 第一类情形，面的某条边是原来直线的子线段 */
					for (edgCod = 0; edgCod < 3; edgCod++)
					{
						if (oriIsZero(iOrtArr[edgCod]) && oriIsZero(jOrtArr[edgCod]))
						{
							/* 确定面片这条边的2个点都在边的内部 */
							iCord = fabs(linep[0][0] - linep[1][0]) > fabs(linep[0][1] - linep[1][1]) ? 0 : 1;
							d1 = linep[0][iCord] > linep[1][iCord] ? linep[0][iCord] : linep[1][iCord];
							d2 = linep[0][iCord] > linep[1][iCord] ? linep[1][iCord] : linep[0][iCord];

							if (facep[edgCod][iCord] >= d2 && facep[edgCod][iCord] <= d1)
							{
								/* 不在延长线上 */
								/* 另外一个点也在线的内部 */
								assert(facep[(edgCod + 1) % 3][iCord] >= d2 && facep[(edgCod + 1) % 3][iCord] <= d1);
								isInt = 2;
								*intCnt = 2;
								intTyp[0] = LTI_INTERSECT_NOD;
								intCod[0] = edgCod;
								intTyp[1] = LTI_INTERSECT_NOD;
								intCod[1] = (edgCod + 1) % 3;

								intPnt[0][0] = facep[edgCod][0];
								intPnt[0][1] = facep[edgCod][1];
								intPnt[1][0] = facep[(edgCod + 1) % 3][0];
								intPnt[1][1] = facep[(edgCod + 1) % 3][1];

								return 1;
							}
							else
							{
								*intTyp = LTI_INTERSECT_NUL;
								*intCnt = 0;
								return 0; /* 不相交 */
							}
						}
					}

					/* 第2类情形，面有1个端点在直线上 （走到这，面不可能有2个端点在直线上）*/
					mOrtArr[0] = ORIENT2D_FUN(linep[0], linep[1], facep[0]);
					mOrtArr[1] = ORIENT2D_FUN(linep[0], linep[1], facep[1]);
					mOrtArr[2] = ORIENT2D_FUN(linep[0], linep[1], facep[2]);
					for (nodCod = 0; nodCod < 3; nodCod++)
					{
						if (oriIsZero(mOrtArr[nodCod]))
						{
							/* 是否在延长线上 */
							iCord = fabs(linep[0][0] - linep[1][0]) > fabs(linep[0][1] - linep[1][1]) ? 0 : 1;
							d1 = linep[0][iCord] > linep[1][iCord] ? linep[0][iCord] : linep[1][iCord];
							d2 = linep[0][iCord] > linep[1][iCord] ? linep[1][iCord] : linep[0][iCord];

							if (facep[nodCod][iCord] < d2 || facep[nodCod][iCord] > d1)
							{
								/* 延长线上，不相交 */
								*intTyp = LTI_INTERSECT_NUL;
								*intCnt = 0;
								return 0;
							}

							/* 不在延长线上，第1个交点即为该顶点 */
							isInt = 1;

							instNodNum = 0;
							twoIntTyp[instNodNum] = LTI_INTERSECT_NOD;
							twoIntCod[instNodNum] = nodCod;
							twoIntPnt[instNodNum][0] = facep[nodCod][0];
							twoIntPnt[instNodNum][1] = facep[nodCod][1];
							instNodNum++;

							/* 看看和nodCod相对的边是否和直线相交 */
							edgCod = (nodCod + 1) % 3;
							if (mOrtArr[edgCod] * mOrtArr[(edgCod + 1) % 3] < 0.0 && iOrtArr[edgCod] * jOrtArr[edgCod] < 0.0)
							{
								twoIntTyp[instNodNum] = LTI_INTERSECT_EDG;
								twoIntCod[instNodNum] = edgCod;

								w1 = fabs(jOrtArr[edgCod]);
								w2 = fabs(iOrtArr[edgCod]);
								w1 = w1 / (w1 + w2);
								w2 = 1.0 - w1;
								twoIntPnt[instNodNum][0] = w1 * linep[0][0] + w2 * linep[1][0];
								twoIntPnt[instNodNum][1] = w1 * linep[0][1] + w2 * linep[1][1];
								instNodNum++;
							}

							break;
						}
					}

					/* 第3类情形，面的2条边和直线可能有交点 */
					if (isInt == 0)
					{
						/* 如果直线的2个端点在面片3条边的某条边不同侧时 */
						for (candEdgNum = 0, edgCod = 0; edgCod < 3; edgCod++)
						{
							/* 注意，此时没有等于0的判断，因为如果等于0，说明一个顶点和这条曲线共线
								* 这和这条边的交点必然也位于在这条边上，这在前面已有判断 */
							if (iOrtArr[edgCod] * jOrtArr[edgCod] < 0.0)
							{
								candEdgCod[candEdgNum++] = edgCod;
							}
						}

						if (candEdgNum > 0)
						{
							for (m = 0, instNodNum = 0; m < candEdgNum; m++)
							{
								edgCod = candEdgCod[m];
								/* 边的两个端点是否异号 */
								if (mOrtArr[edgCod] * mOrtArr[(edgCod + 1) % 3] < 0.0)
								{/* 相交 */
									isInt = 1;
									assert(instNodNum < 2);
									if (instNodNum >= 2)
									{
										//							printf("Error in lin_tri_intersect2d(...). Degenerate triangle.\n");
										*intTyp = LTI_INTERSECT_DEG_FACE;
										*intCnt = 1;
										return 1; /* 返回相交 */
									}
									twoIntTyp[instNodNum] = LTI_INTERSECT_EDG;
									twoIntCod[instNodNum] = edgCod;

									w1 = fabs(jOrtArr[edgCod]);
									w2 = fabs(iOrtArr[edgCod]);
									w1 = w1 / (w1 + w2);
									w2 = 1.0 - w1;
									twoIntPnt[instNodNum][0] = w1 * linep[0][0] + w2 * linep[1][0];
									twoIntPnt[instNodNum][1] = w1 * linep[0][1] + w2 * linep[1][1];
									instNodNum++;
								}
							}
						}
					}

					if (isInt != 0)
					{
						assert(instNodNum == 1 || instNodNum == 2);
						*intCnt = instNodNum;
						for (int iInt = 0; iInt < instNodNum; iInt++) {
							intTyp[iInt] = twoIntTyp[iInt];
							intCod[iInt] = twoIntCod[iInt];
							intPnt[iInt][0] = twoIntPnt[iInt][0];
							intPnt[iInt][1] = twoIntPnt[iInt][1];
						}
					}
					else {
						*intCnt = 0;
					}
				}
			}

			return isInt;
		}

		int pnt_tri_intersect3d_V2(double poinp[3], double facep[3][3], int* intTyp, int* intCod, double intPnt[3], double tolerance)
		{
			auto oriIsZero = [tolerance](double orient) {
				return std::abs(orient) <= tolerance;
				};

			int isInt = 0;
			double q1[2], q2[2], q3[2], q4[2];
			double jOrt, j1Ort, j2Ort, j3Ort;

			/* 注意，我们不允许退化成一条直线的平面存在，因此，如果出现这种情形，都直接返回1 */
			if (proj_four_coplanr_points(facep[0], facep[1], facep[2], poinp, q1, q2, q3, q4) == 1)
			{/* 投影成功 */
				jOrt = orient2d(q1, q2, q3);
				j1Ort = orient2d(q1, q2, q4);
				j2Ort = orient2d(q2, q3, q4);
				j3Ort = orient2d(q3, q1, q4);

				/* 如果j1Ort/j2Ort/j3Ort有2个等于0 (投影成功，则排除3个都等于0的可能)，则交点在点上 */
				if (oriIsZero(j1Ort) && oriIsZero(j2Ort))
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 1;
				}
				else if (oriIsZero(j2Ort) && oriIsZero(j3Ort))
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 2;
				}
				else if (oriIsZero(j3Ort) && oriIsZero(j1Ort))
				{
					isInt = 1;
					*intTyp = LTI_INTERSECT_NOD;
					*intCod = 0;
				}

				if (isInt == 0)
				{
					/* 如果j1Ort/j2Ort/j3Ort只有一个等于0，且另外两个与jOrt同号，则交点在边上 */
					if (oriIsZero(j1Ort) && j2Ort * jOrt > 0.0 && j3Ort * jOrt > 0.0)
					{
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 0;
					}
					if (oriIsZero(j2Ort) && j3Ort * jOrt > 0.0 && j1Ort * jOrt > 0.0)
					{
						/* 走到这一步，j3Ort和j1Ort的值必然和jOrt同号 */
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 1;
					}
					if (oriIsZero(j3Ort) && j1Ort * jOrt > 0.0 && j2Ort * jOrt > 0.0)
					{
						/* 走到这一步，j1Ort和j2Ort的值必然和jOrt同号 */
						isInt = 1;
						*intTyp = LTI_INTERSECT_EDG;
						*intCod = 2;
					}
				}

				if (isInt == 0)
				{
					/* 如果j1Ort/j2Ort/j3Ort都与jOrt同号，则交点在面内 */
					if (j1Ort * jOrt > 0.0 && j2Ort * jOrt > 0.0 && j3Ort * jOrt > 0.0)
					{
						isInt = 1;
						*intTyp = LTI_INTERSECT_FAC;
					}
				}

				/* 如果相交，保存交点 */
				if (isInt && intPnt)
				{
					intPnt[0] = poinp[0];
					intPnt[1] = poinp[1];
					intPnt[2] = poinp[2];
				}
			}

			return isInt;
		}
	}
}