package fq.geom.nurbs;

/**
 * This class contains Algorithms from <i>The NURBS Book</i>.
 * 
 * @author FuQiang
 * 
 */
public class Algorithms {
	/**
	 * Determine the knot span index
	 * 
	 * @param u -
	 *            parameter value
	 * @param deg -
	 *            degree of B-Spline basis function
	 * @param knots -
	 *            nonperiodic knot vector
	 * 
	 * @return the kont span index
	 */
	static public int FindSpan(double u, int deg, double[] knots) {
		int low, high, mid;

		int n = knots.length - deg - 2;
		// special case
		if (u >= knots[n + 1])
			return n;
		if (u <= knots[deg])
			return deg;

		// do binary search
		low = deg;
		high = n + 1;
		mid = (low + high) / 2;
		while (u < knots[mid] || u >= knots[mid + 1]) {
			if (u < knots[mid])
				high = mid;
			else
				low = mid;
			mid = (low + high) / 2;
		}
		return (mid);
	}

	/**
	 * @param u -
	 *            parameter value
	 * @param span -
	 *            kont span index
	 * @param deg -
	 *            degree of B-Spline basis function
	 * @param knots -
	 *            nonperiodic knot vector
	 * 
	 * @return all the nonvanishing basis functions: N[0],...,N[p].
	 */
	static public double[] BasisFuns(double u, int span, int deg, double[] knots) {
		double[] N = new double[deg + 1];
		N[0] = 1.0;

		double[] left = new double[deg + 1];
		double[] right = new double[deg + 1];
		double saved, temp;
		for (int j = 1; j <= deg; j++) {
			left[j] = u - knots[span + 1 - j];
			right[j] = knots[span + j] - u;
			saved = 0.0;
			for (int r = 0; r < j; r++) {
				temp = N[r] / (right[r + 1] + left[j - r]);
				N[r] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			N[j] = saved;
		}
		return N;
	}

	/**
	 * Generate nonperiodic knot vectors, within range [0,1]
	 * 
	 * @param n -
	 *            (n+1) is the number of control points
	 * @param deg -
	 *            degree of the B-spline basis functions
	 * @return Generated nonperiodic knot vectors
	 */
	static public double[] GenKnots(int n, int deg) {
		int m = n + deg + 1;
		double inv = 1 / (double) (n - deg + 1);

		double[] U = new double[m + 1];
		for (int i = 0; i <= m; i++) {
			if (i <= deg)
				U[i] = 0.0;
			else if (i > deg && i <= n)
				U[i] = (double) (i - deg) * inv;
			else
				U[i] = 1.0;
		}
		return U;
	}

	/**
	 * Generate nonperiodic knot vectors, within range [0,1]
	 * 
	 * @param n -
	 *            (n+1) is the number of control points
	 * @param deg -
	 *            degree of the B-spline basis functions
	 * @return Generated nonperiodic knot vectors
	 */
	static public boolean IsKnotsValid(int n, int deg, double[] newKnots) {
		// m = n + p + 1
		int m = n + deg + 1;
		
		if(newKnots.length != m+1)
			return false;
		
		if(newKnots[m] <= newKnots[0])
			return false;
			
		for (int i = 0; i < m; i++) {
			if(newKnots[i+1] < newKnots[i])
				return false;
		}
		//for (int i =0; i<deg;i++) {
		//	if(newKnots[i] != newKnots[i+1])
		//		return false;
		//	if(newKnots[m - i] != newKnots[m-i-1])
		//		return false;
		//}
		return true;
	}

	/**
	 * Compute the nonzero basis functions and their derivatives, up to and
	 * including the nth derivative (n<=deg).
	 * 
	 * @param u -
	 *            parameter value
	 * @param n -
	 *            nth derivative
	 * @param span -
	 *            knot vector span index
	 * @param deg -
	 *            degree of basis function
	 * @param U -
	 *            knot vector
	 * @return A two-dimensional array, Ders[][]<br>
	 *         Ders[k][j] is the kth derivative of the function
	 *         N(span-deg+j,deg),<br>
	 *         where k:[0,n] and j:[0,p]
	 */
	static public double[][] DersBasisFuns(double u, int n, int span, int deg,
			double[] U) {
		double[][] ders = new double[n + 1][deg + 1];

		double[] left = new double[deg + 1];
		double[] right = new double[deg + 1];
		double[][] ndu = new double[deg + 1][deg + 1];
		double saved, temp;
		int j, r;

		ndu[0][0] = 1.0;
		for (j = 1; j <= deg; j++) {
			left[j] = u - U[span + 1 - j];
			right[j] = U[span + j] - u;
			saved = 0.0;
			for (r = 0; r < j; r++) {
				// Lower triangle
				ndu[j][r] = right[r + 1] + left[j - r];
				temp = ndu[r][j - 1] / ndu[j][r];
				// Upper triangle
				ndu[r][j] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			ndu[j][j] = saved;
		}

		// load the basis functions
		for (j = 0; j <= deg; j++) {
			ders[0][j] = ndu[j][deg];
		}

		// compute the derivatives, Eq.[2.9] the NURBS Book
		double[][] a = new double[2][deg + 1];
		for (r = 0; r <= deg; r++) { // loop over function index
			int s1, s2;
			s1 = 0;
			s2 = 1;
			a[0][0] = 1.0;
			// compute the kth derivative
			for (int k = 1; k <= n; k++) {
				double d;
				int rk, pk, j1, j2;

				d = 0.0;
				rk = r - k;
				pk = deg - k;
				if (r >= k) {
					a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
					d = a[s2][0] * ndu[rk][pk];
				}

				if (rk >= -1) {
					j1 = 1;
				} else {
					j1 = -rk;
				}
				if (r - 1 <= pk) {
					j2 = k - 1;
				} else {
					j2 = deg - r;
				}

				for (j = j1; j <= j2; j++) {
					a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
					d += a[s2][j] * ndu[rk + j][pk];
				}
				if (r <= pk) {
					a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
					d += a[s2][j] * ndu[r][pk];
				}
				ders[k][r] = d;
				j = s1;
				s1 = s2;
				s2 = j; // switch rows
			}
		}
		// Multiply through by the correct factors
		// Eq [2.9] of the NURBS Book
		r = deg;
		for (int k = 1; k <= n; k++) {
			for (j = 0; j <= deg; j++)
				ders[k][j] *= r;
			r *= (deg - k);
		}
		return ders;
	}

	/**
	 * Compute binomial coefficients, and store in an array
	 * 
	 * @return An array of binomial coefficients, Bin[K][I]
	 * 
	 * <pre>
	 *                     k!
	 *    Bin[k][i] =  -----------
	 *                  i! (k-i)!
	 * </pre>
	 */
	public static double[][] BinomialCoef(int K, int I) {
		double Bin[][] = new double[K + 1][I + 1];
		int k, i;
		// Setup the first line
		Bin[0][0] = 1.0;
		for (i = I; i > 0; --i)
			Bin[0][i] = 0.0;
		// Setup the other lines
		for (k = 0; k < K; k++) {
			Bin[k + 1][0] = 1.0;
			for (i = 1; i <= I; i++) {
				if (k + 1 < i)
					Bin[k][i] = 0.0;
				else
					Bin[k + 1][i] = Bin[k][i] + Bin[k][i - 1];
			}
		}
		return Bin;
	}
}
