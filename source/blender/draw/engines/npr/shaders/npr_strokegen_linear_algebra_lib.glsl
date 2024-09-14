#if !defined(INCLUDE_NPR_STROKEGEN_LINEAR_ALGEBRA_LIB_H)
#define INCLUDE_NPR_STROKEGEN_LINEAR_ALGEBRA_LIB_H



/* 3x3 Eigen Solver from https://github.com/hjwdzh/Fluid3D/blob/master/FluidCS11.hlsl 
 * ---------------------------------------------------------------------------------------- */

// Constants
#define M_SQRT3    1.73205081f   // sqrt(3)
#define FLT_EPSILON     1.192092896e-07f

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 

float atan2_hlsl(float x, float y)
{ // https://community.khronos.org/t/hlsl-to-glsl-atan2-function/53881/4
    return atan(x, y);
}

// ----------------------------------------------------------------------------
void dsyevc3(mat3 A, inout vec3 w)
	// ----------------------------------------------------------------------------
	// Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
	// analytical algorithm.
	// Only the diagonal and upper triangular parts of A are accessed. The access
	// is read-only.
	// ----------------------------------------------------------------------------
	// Parameters:
	//   A: The symmetric input matrix
	//   w: Storage buffer for eigenvalues
	// ----------------------------------------------------------------------------
	// Return value:
	//   0: Success
	//  -1: Error
	// ----------------------------------------------------------------------------
{
	float m, c1, c0;

	// Determine coefficients of characteristic poynomial. We write
	//       | a   d   f  |
	//  A =  | d*  b   e  |
	//       | f*  e*  c  |
	float de = A[0][1] * A[1][2];                                    // d * e
	float dd = SQR(A[0][1]);                                         // d^2
	float ee = SQR(A[1][2]);                                         // e^2
	float ff = SQR(A[0][2]);                                         // f^2
	m  = A[0][0] + A[1][1] + A[2][2];
	c1 = (A[0][0]*A[1][1] + A[0][0]*A[2][2] + A[1][1]*A[2][2])        // a*b + a*c + b*c - d^2 - e^2 - f^2
		- (dd + ee + ff);
	c0 = A[2][2]*dd + A[0][0]*ee + A[1][1]*ff - A[0][0]*A[1][1]*A[2][2]
	- 2.0f * A[0][2]*de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

	float p, sqrt_p, q, c, s, phi;
	p = SQR(m) - 3.0f*c1;
	q = m*(p - (3.0f/2.0f)*c1) - (27.0f/2.0f)*c0;
	sqrt_p = sqrt(abs(p));

	phi = 27.0f * ( 0.25f*SQR(c1)*(p - c1) + c0*(q + 27.0f/4.0f*c0));
	phi = (1.0f/3.0f) * atan2_hlsl(sqrt(abs(phi)), q);

	c = sqrt_p*cos(phi);
	s = (1.0f/M_SQRT3)*sqrt_p*sin(phi);

	w[1]  = (1.0f/3.0f)*(m - c);
	w[2]  = w[1] + s;
	w[0]  = w[1] + c;
	w[1] -= s;
}

// ----------------------------------------------------------------------------
void dsyevv3(mat3 A, inout mat3 Q, inout vec3 w)
	// ----------------------------------------------------------------------------
	// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
	// matrix A using Cardano's method for the eigenvalues and an analytical
	// method based on vector cross products for the eigenvectors.
	// Only the diagonal and upper triangular parts of A need to contain meaningful
	// values. However, all of A may be used as temporary storage and may hence be
	// destroyed.
	// ----------------------------------------------------------------------------
	// Parameters:
	//   A: The symmetric input matrix
	//   Q: Storage buffer for eigenvectors
	//   w: Storage buffer for eigenvalues
	// ----------------------------------------------------------------------------
	// Return value:
	//   0: Success
	//  -1: Error
	// ----------------------------------------------------------------------------
	// Dependencies:
	//   dsyevc3()
	// ----------------------------------------------------------------------------
	// Version history:
	//   v1.1 (12 Mar 2012): Removed access to lower triangualr part of A
	//     (according to the documentation, only the upper triangular part needs
	//     to be filled)
	//   v1.0f: First released version
	// ----------------------------------------------------------------------------
{
#ifndef EVALS_ONLY
	float norm;          // Squared norm or inverse norm of current eigenvector
	float n0, n1;        // Norm of first and second columns of A
	float n0tmp, n1tmp;  // "Templates" for the calculation of n0/n1 - saves a few FLOPS
	float thresh;        // Small number used as threshold for floating point comparisons
	float error;         // Estimated maximum roundoff error in some steps
	float wmax;          // The eigenvalue of maximum modulus
	float f, t;          // Intermediate storage
	int i, j;             // Loop counters
#endif

	// Calculate eigenvalues
	dsyevc3(A, w);

#ifndef EVALS_ONLY
	wmax = abs(w[0]);
	if ((t=abs(w[1])) > wmax)
		wmax = t;
	if ((t=abs(w[2])) > wmax)
		wmax = t;
	thresh = SQR(8.0f * FLT_EPSILON * wmax);

	// Prepare calculation of eigenvectors
	n0tmp   = SQR(A[0][1]) + SQR(A[0][2]);
	n1tmp   = SQR(A[0][1]) + SQR(A[1][2]);
	Q[0][1] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
	Q[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0];
	Q[2][1] = SQR(A[0][1]);

	// Calculate first eigenvector by the formula
	//   v[0] = (A - w[0]).e1 x (A - w[0]).e2
	A[0][0] -= w[0];
	A[1][1] -= w[0];
	Q[0][0] = Q[0][1] + A[0][2]*w[0];
	Q[1][0] = Q[1][1] + A[1][2]*w[0];
	Q[2][0] = A[0][0]*A[1][1] - Q[2][1];
	norm    = SQR(Q[0][0]) + SQR(Q[1][0]) + SQR(Q[2][0]);
	n0      = n0tmp + SQR(A[0][0]);
	n1      = n1tmp + SQR(A[1][1]);
	error   = n0 * n1;

	if (n0 <= thresh)         // If the first column is zero, then (1,0,0) is an eigenvector
	{
		Q[0][0] = 1.0f;
		Q[1][0] = 0.0f;
		Q[2][0] = 0.0f;
	}
	else if (n1 <= thresh)    // If the second column is zero, then (0,1,0) is an eigenvector
	{
		Q[0][0] = 0.0f;
		Q[1][0] = 1.0f;
		Q[2][0] = 0.0f;
	}
	else if (norm < SQR(64.0f * FLT_EPSILON) * error)
	{                         // If angle between A[0] and A[1] is too small, don't use
		t = SQR(A[0][1]);       // cross product, but calculate v ~ (1, -A0/A1, 0)
		f = -A[0][0] / A[0][1];
		if (SQR(A[1][1]) > t)
		{
			t = SQR(A[1][1]);
			f = -A[0][1] / A[1][1];
		}
		if (SQR(A[1][2]) > t)
			f = -A[0][2] / A[1][2];
		norm    = 1.0f/sqrt(1 + SQR(f));
		Q[0][0] = norm;
		Q[1][0] = f * norm;
		Q[2][0] = 0.0f;
	}
	else                      // This is the standard branch
	{
		norm = sqrt(1.0f / norm);
		for (j=0; j < 3; j++)
			Q[j][0] = Q[j][0] * norm;
	}


	// Prepare calculation of second eigenvector
	t = w[0] - w[1];
	if (abs(t) > 8.0f * FLT_EPSILON * wmax)
	{
		// For non-degenerate eigenvalue, calculate second eigenvector by the formula
		//   v[1] = (A - w[1]).e1 x (A - w[1]).e2
		A[0][0] += t;
		A[1][1] += t;
		Q[0][1]  = Q[0][1] + A[0][2]*w[1];
		Q[1][1]  = Q[1][1] + A[1][2]*w[1];
		Q[2][1]  = A[0][0]*A[1][1] - Q[2][1];
		norm     = SQR(Q[0][1]) + SQR(Q[1][1]) + SQR(Q[2][1]);
		n0       = n0tmp + SQR(A[0][0]);
		n1       = n1tmp + SQR(A[1][1]);
		error    = n0 * n1;

		if (n0 <= thresh)       // If the first column is zero, then (1,0,0) is an eigenvector
		{
			Q[0][1] = 1.0f;
			Q[1][1] = 0.0f;
			Q[2][1] = 0.0f;
		}
		else if (n1 <= thresh)  // If the second column is zero, then (0,1,0) is an eigenvector
		{
			Q[0][1] = 0.0f;
			Q[1][1] = 1.0f;
			Q[2][1] = 0.0f;
		}
		else if (norm < SQR(64.0f * FLT_EPSILON) * error)
		{                       // If angle between A[0] and A[1] is too small, don't use
			t = SQR(A[0][1]);     // cross product, but calculate v ~ (1, -A0/A1, 0)
			f = -A[0][0] / A[0][1];
			if (SQR(A[1][1]) > t)
			{
				t = SQR(A[1][1]);
				f = -A[0][1] / A[1][1];
			}
			if (SQR(A[1][2]) > t)
				f = -A[0][2] / A[1][2];
			norm    = 1.0f/sqrt(1 + SQR(f));
			Q[0][1] = norm;
			Q[1][1] = f * norm;
			Q[2][1] = 0.0f;
		}
		else
		{
			norm = sqrt(1.0f / norm);
			for (j=0; j < 3; j++)
				Q[j][1] = Q[j][1] * norm;
		}
	}
	else
	{
		// For degenerate eigenvalue, calculate second eigenvector according to
		//   v[1] = v[0] x (A - w[1]).e[i]
		//   
		// This would really get to complicated if we could not assume all of A to
		// contain meaningful values.
		A[1][0]  = A[0][1];
		A[2][0]  = A[0][2];
		A[2][1]  = A[1][2];
		A[0][0] += w[0];
		A[1][1] += w[0];

		for (i=0; i < 3; i++)
		{
			A[i][i] -= w[1];
			n0       = SQR(A[0][i]) + SQR(A[1][i]) + SQR(A[2][i]);
			if (n0 > thresh)
			{
				Q[0][1]  = Q[1][0]*A[2][i] - Q[2][0]*A[1][i];
				Q[1][1]  = Q[2][0]*A[0][i] - Q[0][0]*A[2][i];
				Q[2][1]  = Q[0][0]*A[1][i] - Q[1][0]*A[0][i];
				norm     = SQR(Q[0][1]) + SQR(Q[1][1]) + SQR(Q[2][1]);
				if (norm > SQR(256.0f * FLT_EPSILON) * n0) // Accept cross product only if the angle between
				{                                         // the two vectors was not too small
					norm = sqrt(1.0f / norm);
					for (j=0; j < 3; j++)
						Q[j][1] = Q[j][1] * norm;
					break;
				}
			}
		}

		if (i == 3)    // This means that any vector orthogonal to v[0] is an EV.
		{
			for (j=0; j < 3; j++)
				if (Q[j][0] != 0.0f)                                   // Find nonzero element of v[0] ...
				{                                                     // ... and swap it with the next one
					norm          = 1.0f / sqrt(SQR(Q[j][0]) + SQR(Q[(j+1)%3][0]));
					Q[j][1]       = Q[(j+1)%3][0] * norm;
					Q[(j+1)%3][1] = -Q[j][0] * norm;
					Q[(j+2)%3][1] = 0.0f;
					break;
				}
		}
	}


	// Calculate third eigenvector according to
	//   v[2] = v[0] x v[1]
	Q[0][2] = Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1];
	Q[1][2] = Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1];
	Q[2][2] = Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1];
#endif
}










// -----------------------------------------------------------------------------------
// "A Robust Eigensolver for 3 × 3 Symmetric Matrices"
// https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf#page=10.55

#define T float // save my time for translating the "T"s to "float"s

void ComputeOrthogonalComplement(vec3 W, inout vec3 U, inout vec3 V) 
{
	// Robustly compute a right-handed orthonormal set { U, V, W }.
	// The vector W is guaranteed to be unit-length, in which case
	// there is no need to worry about a division by zero when
	// computing invLength.
	T invLength;
	if (abs(W[0]) > abs(W[1]))
	{
		// The component of maximum absolute value is either W[0]
		// or W[2].
		invLength = 1.0f / sqrt(W[0] * W[0] + W[2] * W[2]);
		U = vec3( -W[2] * invLength, .0f, +W[0] * invLength );
	}
	else
	{
		// The component of maximum absolute value is either W[1]
		// or W[2].
		invLength = 1.0f / sqrt(W[1] * W[1] + W[2] * W[2]);
		U = vec3( .0f, +W[2] * invLength, -W[1] * invLength );
	}
	V = cross(W, U);
}

void ComputeEigenvector0(
	T a00, T a01, T a02, T a11, T a12, T a22, T eval0, 
	inout vec3 evec0
) {
	// Compute a unit-length eigenvector for eigenvalue[i0].  The
	// matrix is rank 2, so two of the rows are linearly independent.
	// For a robust computation of the eigenvector, select the two
	// rows whose cross product has largest length of all pairs of
	// rows.
	vec3 row0 = vec3( a00 - eval0, a01, a02 );
	vec3 row1 = vec3( a01, a11 - eval0, a12 );
	vec3 row2 = vec3( a02, a12, a22 - eval0 );
	vec3 r0xr1 = cross(row0, row1);
	vec3 r0xr2 = cross(row0, row2);
	vec3 r1xr2 = cross(row1, row2);
	T d0 = dot(r0xr1, r0xr1);
	T d1 = dot(r0xr2, r0xr2);
	T d2 = dot(r1xr2, r1xr2);

	T dmax = d0;
	int imax = 0;
	if (d1 > dmax)
	{
		dmax = d1;
		imax = 1;
	}
	if (d2 > dmax)
	{
		imax = 2;
	}

	if (imax == 0)
	{
		evec0 = r0xr1 / sqrt(d0);
	}
	else if (imax == 1)
	{
		evec0 = r0xr2 / sqrt(d1);
	}
	else
	{
		evec0 = r1xr2 / sqrt(d2);
	}
}

void ComputeEigenvector1(
	T a00, T a01, T a02, T a11, T a12, T a22,
	vec3 evec0, T eval1, 
	inout vec3 evec1
) {
	// Robustly compute a right-handed orthonormal set
	// { U, V, evec0 }.
	vec3 U, V;
	ComputeOrthogonalComplement(evec0, U, V);

	// Let e be eval1 and let E be a corresponding eigenvector which
	// is a solution to the linear system (A - e*I)*E = 0.  The matrix
	// (A - e*I) is 3x3, not invertible (so infinitely many
	// solutions), and has rank 2 when eval1 and eval are different.
	// It has rank 1 when eval1 and eval2 are equal.  Numerically, it
	// is difficult to compute robustly the rank of a matrix.  Instead,
	// the 3x3 linear system is reduced to a 2x2 system as follows.
	// Define the 3x2 matrix J = [U V] whose columns are the U and V
	// computed previously.  Define the 2x1 vector X = J*E.  The 2x2
	// system is 0 = M * X = (J^T * (A - e*I) * J) * X where J^T is
	// the transpose of J and M = J^T * (A - e*I) * J is a 2x2 matrix.
	// The system may be written as
	//     +-                        -++-  -+       +-  -+
	//     | U^T*A*U - e  U^T*A*V     || x0 | = e * | x0 |
	//     | V^T*A*U      V^T*A*V - e || x1 |       | x1 |
	//     +-                        -++   -+       +-  -+
	// where X has row entries x0 and x1.

	vec3 AU =
	{
		a00 * U[0] + a01 * U[1] + a02 * U[2],
		a01 * U[0] + a11 * U[1] + a12 * U[2],
		a02 * U[0] + a12 * U[1] + a22 * U[2]
	};

	vec3 AV =
	{
		a00 * V[0] + a01 * V[1] + a02 * V[2],
		a01 * V[0] + a11 * V[1] + a12 * V[2],
		a02 * V[0] + a12 * V[1] + a22 * V[2]
	};

	T m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eval1;
	T m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
	T m11 = V[0] * AV[0] + V[1] * AV[1] + V[2] * AV[2] - eval1;

	// For robustness, choose the largest-length row of M to compute
	// the eigenvector.  The 2-tuple of coefficients of U and V in the
	// assignments to eigenvector[1] lies on a circle, and U and V are
	// unit length and perpendicular, so eigenvector[1] is unit length
	// (within numerical tolerance).
	T absM00 = abs(m00);
	T absM01 = abs(m01);
	T absM11 = abs(m11);
	T maxAbsComp;
	if (absM00 >= absM11)
	{
		maxAbsComp = max(absM00, absM01);
		if (maxAbsComp > .0f)
		{
			if (absM00 >= absM01)
			{
				m01 /= m00;
				m00 = 1.0f / sqrt(1.0f + m01 * m01);
				m01 *= m00;
			}
			else
			{
				m00 /= m01;
				m01 = 1.0f / sqrt(1.0f + m00 * m00);
				m00 *= m01;
			}
			evec1 = m01 * U - m00 * V;
		}
		else
		{
			evec1 = U;
		}
	}
	else
	{
		maxAbsComp = max(absM11, absM01);
		if (maxAbsComp > .0f)
		{
			if (absM11 >= absM01)
			{
				m01 /= m11;
				m11 = 1.0f / sqrt(1.0f + m01 * m01);
				m01 *= m11;
			}
			else
			{
				m11 /= m01;
				m01 = 1.0f / sqrt(1.0f + m11 * m11);
				m11 *= m01;
			}
			evec1 = m11 * U - m01 * V;
		}
		else
		{
			evec1 = U;
		}
	}
}

void NISymmetricEigensolver3x3_solve(
	T a00, T a01, T a02, T a11, T a12, T a22,
	// uint sortType, 
	inout vec3 eval, inout mat3 evec
) 
{
	// Precondition the matrix by factoring out the maximum absolute
	// value of the components.  This guards against floating-point
	// overflow when computing the eigenvalues.
	T max0 = max(abs(a00), abs(a01));
	T max1 = max(abs(a02), abs(a11));
	T max2 = max(abs(a12), abs(a22));
	T maxAbsElement = max(max(max0, max1), max2);
	if (maxAbsElement == .0f)
	{
		// A is the zero matrix.
		eval[0] = .0f;
		eval[1] = .0f;
		eval[2] = .0f;
		evec[0] = vec3( 1.0f, .0f, .0f );
		evec[1] = vec3( .0f, 1.0f, .0f );
		evec[2] = vec3( .0f, .0f, 1.0f );
		return;
	}

	T invMaxAbsElement = 1.0f / maxAbsElement;
	a00 *= invMaxAbsElement;
	a01 *= invMaxAbsElement;
	a02 *= invMaxAbsElement;
	a11 *= invMaxAbsElement;
	a12 *= invMaxAbsElement;
	a22 *= invMaxAbsElement;

	T norm = a01 * a01 + a02 * a02 + a12 * a12;
	if (norm > T(0))
	{
		// Compute the eigenvalues of A.

		// In the PDF mentioned previously, B = (A - q*I)/p, where
		// q = tr(A)/3 with tr(A) the trace of A (sum of the diagonal
		// entries of A) and where p = sqrt(tr((A - q*I)^2)/6).
		T q = (a00 + a11 + a22) / 3.0f;

		// The matrix A - q*I is represented by the following, where
		// b00, b11 and b22 are computed after these comments,
		//   +-           -+
		//   | b00 a01 a02 |
		//   | a01 b11 a12 |
		//   | a02 a12 b22 |
		//   +-           -+
		T b00 = a00 - q;
		T b11 = a11 - q;
		T b22 = a22 - q;

		// The is the variable p mentioned in the PDF.
		T p = sqrt((b00 * b00 + b11 * b11 + b22 * b22 + norm * 2.0f) / 6.0f);

		// We need det(B) = det((A - q*I)/p) = det(A - q*I)/p^3.  The
		// value det(A - q*I) is computed using a cofactor expansion
		// by the first row of A - q*I.  The cofactors are c00, c01
		// and c02 and the determinant is b00*c00 - a01*c01 + a02*c02.
		// The det(B) is then computed finally by the division
		// with p^3.
		T c00 = b11 * b22 - a12 * a12;
		T c01 = a01 * b22 - a12 * a02;
		T c02 = a01 * a12 - b11 * a02;
		T det = (b00 * c00 - a01 * c01 + a02 * c02) / (p * p * p);

		// The halfDet value is cos(3*theta) mentioned in the PDF. The
		// acos(z) function requires |z| <= 1, but will fail silently
		// and return NaN if the input is larger than 1 in magnitude.
		// To avoid this problem due to rounding errors, the halfDet
		// value is clamped to [-1,1].
		T halfDet = det * 0.5f;
		halfDet = min(max(halfDet, -1.0f), 1.0f);

		// The eigenvalues of B are ordered as
		// beta0 <= beta1 <= beta2.  The number of digits in
		// twoThirdsPi is chosen so that, whether float or double,
		// the floating-point number is the closest to theoretical
		// 2*pi/3.
		T angle = acos(halfDet) / 3.0f;
		T twoThirdsPi = 2.09439510239319549f;
		T beta2 = cos(angle) * 2.0f;
		T beta0 = cos(angle + twoThirdsPi) * 2.0f;
		T beta1 = -(beta0 + beta2);

		// The eigenvalues of A are ordered as
		// alpha0 <= alpha1 <= alpha2.
		eval[0] = q + p * beta0;
		eval[1] = q + p * beta1;
		eval[2] = q + p * beta2;

		// Compute the eigenvectors so that the set
		// {evec[0], evec[1], evec[2]} is right handed and
		// orthonormal.
		if (halfDet >= .0f)
		{
			ComputeEigenvector0(a00, a01, a02, a11, a12, a22, eval[2], evec[2]);
			ComputeEigenvector1(a00, a01, a02, a11, a12, a22, evec[2], eval[1], evec[1]);
			evec[0] = cross(evec[1], evec[2]);
		}
		else
		{
			ComputeEigenvector0(a00, a01, a02, a11, a12, a22, eval[0], evec[0]);
			ComputeEigenvector1(a00, a01, a02, a11, a12, a22, evec[0], eval[1], evec[1]);
			evec[2] = cross(evec[0], evec[1]);
		}
	}
	else
	{
		// The matrix is diagonal.
		eval[0] = a00;
		eval[1] = a11;
		eval[2] = a22;
		evec[0] = vec3( 1.0f, .0f, .0f );
		evec[1] = vec3( .0f, 1.0f, .0f );
		evec[2] = vec3( .0f, .0f, 1.0f );
	}

	// The preconditioning scaled the matrix A, which scales the
	// eigenvalues.  Revert the scaling.
	eval[0] *= maxAbsElement;
	eval[1] *= maxAbsElement;
	eval[2] *= maxAbsElement;

	// SortEigenstuff<T>()(sortType, true, eval, evec);
}
#undef T


#endif 


