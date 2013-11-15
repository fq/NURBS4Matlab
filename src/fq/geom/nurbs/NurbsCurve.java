package fq.geom.nurbs;

import fq.geom.*;

/**
 * NURBS curve class<br>
 * This class is used to represent and manipulate NURBS curve.
 * 
 * @author FuQiang
 * 
 */
public class NurbsCurve extends ParametricCurve {
	double[] knots; // knot vector

	Point3DH[] controlPoints; // control points

	int degree; // degree

	/**
	 * Initialize the NURRBS curve
	 * 
	 * @param cp -
	 *            Array of N control points<br>
	 *            cp(i,1) -> x coordinates of the ith control points<br>
	 *            cp(i,2) -> y coordinates of the ith control points<br>
	 *            cp(i,3) -> z coordinates of the ith control points(optional,
	 *            default is 0.0)<br>
	 *            cp(i,4) -> weights of ith control points(optional, default is
	 *            1.0)<br>
	 *            where, 0 <= i <= N;
	 * @param deg -
	 *            Degree of the NURBS curve
	 */
	public NurbsCurve(double[][] cp, int deg) {
		init(cp, deg);
	}

	public NurbsCurve(double[][] cp) {
		init(cp, 3);
	}

	/**
	 * Initialize the NURBS curve
	 */
	void init(double[][] cp, int deg) {
		if (deg < 1)
			throw new GeomException("Degree must >= 1!");
		int n = cp.length - 1;
		if (n < 1)
			throw new GeomException("at least two control points!");
		if (deg >= n)
			deg = n;
		controlPoints = new Point3DH[n + 1];
		for (int i = 0; i <= n; i++) {
			controlPoints[i] = new Point3DH(cp[i]);
		}
		degree = deg;
		knots = Algorithms.GenKnots(n, degree);
	}

	/**
	 * Determine the knot span index at parameter u
	 * 
	 * @param u -
	 *            the parametric value
	 * @return the knot span index (zero based)
	 */
	public int findSpan(double u) {
		return Algorithms.FindSpan(u, degree, knots);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see fq.geom.ParametricCurve#getPointAt(double)
	 */
	public Point3D getPointAt(double u) {
		if (u < minKnot() || u > maxKnot())
			throw new GeomException("Input parameter out of range!");
		int span = findSpan(u);
		double[] N = Algorithms.BasisFuns(u, span, degree, knots);
		double xw = 0.0, yw = 0.0, zw = 0.0, ww = 0.0;
		for (int i = 0; i <= degree; i++) {
			xw += N[i] * controlPoints[span - degree + i].XW();
			yw += N[i] * controlPoints[span - degree + i].YW();
			zw += N[i] * controlPoints[span - degree + i].ZW();
			ww += N[i] * controlPoints[span - degree + i].W();
		}
		return new Point3D(xw / ww, yw / ww, zw / ww);
	}

	/**
	 * Compuate derivates on the NURBS curve at homogeneous coordinates
	 */
	double[][] getDerivsAtH(double u, int d) {
		if (u < minKnot() || u > maxKnot())
			throw new GeomException("Input parameter out of range!");
		double[][] dersH = new double[d + 1][4];
		int du = Math.min(degree, d);
		int k;

		int span = findSpan(u);
		double[][] nders = Algorithms.DersBasisFuns(u, du, span, degree, knots);
		for (k = 0; k <= du; k++) {
			for (int j = 0; j <= degree; j++) {
				dersH[k][0] += nders[k][j]
						* controlPoints[span - degree + j].XW();
				dersH[k][1] += nders[k][j]
						* controlPoints[span - degree + j].YW();
				dersH[k][2] += nders[k][j]
						* controlPoints[span - degree + j].ZW();
				dersH[k][3] += nders[k][j]
						* controlPoints[span - degree + j].W();
			}
		}
		return dersH;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see fq.geom.ParametricCurve#getDerivsAt(double, int)
	 */
	public Vector3D[] getDerivsAt(double u, int d) {
		if (u < minKnot() || u > maxKnot())
			throw new GeomException("Input parameter out of range!");

		int du = Math.min(degree, d);
		double[][] dH = getDerivsAtH(u, du);
		int k, i;
		// double[][] CK = new double[d + 1][3];
		Vector3D[] CK = new Vector3D[d + 1];
		double[] v = new double[3];
		double[][] Bin = Algorithms.BinomialCoef(degree, degree);

		for (k = 0; k <= du; k++) {
			v[0] = dH[k][0];
			v[1] = dH[k][1];
			v[2] = dH[k][2];
			for (i = 1; i <= k; i++) {
				v[0] = v[0] - (Bin[k][i] * dH[i][3]) * CK[k - i].X();
				v[1] = v[1] - (Bin[k][i] * dH[i][3]) * CK[k - i].Y();
				v[2] = v[2] - (Bin[k][i] * dH[i][3]) * CK[k - i].Z();
			}
			CK[k] = new Vector3D(v[0] / dH[0][3], v[1] / dH[0][3], v[2]
					/ dH[0][3]);
			// CK[k][0] = v[0] / dH[0][3];
			// CK[k][1] = v[1] / dH[0][3];
			// CK[k][2] = v[2] / dH[0][3];
		}
		return CK;
	}

	/**
	 * Returns the minimum value of the knot vector
	 * 
	 * @return the minimum value of the knot vector
	 */
	public double minKnot() {
		return knots[0];
	}

	/**
	 * Returns the maximum value of the knot vector
	 * 
	 * @return the maximum value of the knot vector
	 */
	public double maxKnot() {
		return knots[knots.length - 1];
	}

	/**
	 * Returns the degree of the NURBS curve
	 */
	public int getDegree() {
		return degree;
	}

	// /**
	// * Sets the degree of the NURBS curve
	// *
	// * @param deg -
	// * New degree for the curve
	// */
	// public void setDegree(int deg) {
	// if (deg < 1)
	// throw new GeomException("Degree must >= 1!");
	//
	// int n = controlPoints.length - 1;
	// if (deg >= n)
	// deg = n;
	// this.degree = deg;
	// knots = Algorithms.GenKnots(n, degree);
	// }

	/**
	 * Returns the number of Control Points
	 * 
	 * @return The number of control points
	 */
	public int getNumberOfControlPoints() {
		return controlPoints.length;
	}

	// /**
	// * Returns the Control Point at given index.<br>
	// * (This method is designed for using in MATLAB)
	// *
	// * @param index -
	// * The index of the point to retrieve
	// */
	// public double[] getControlPoint(int index) {
	// if (index < 0) {
	// index = 0;
	// }
	// if (index > controlPoints.length) {
	// index = controlPoints.length;
	// }
	// return controlPoints[index].toArray();
	// }

	/**
	 * Returns the control points of the NURBS curve
	 * 
	 * @return An (Nx4) array CP of control points<br>
	 *         &nbspN - Number of control points.<br>
	 *         &nbsp4 - Coordinates of control point, including weights.<br>
	 * 
	 */
	public double[][] getControlPoints() {
		double[][] p = new double[controlPoints.length][];
		for (int i = 0; i < controlPoints.length; i++)
			p[i] = controlPoints[i].toArray();
		return p;
	}

	/**
	 * Returns the knot vector
	 * 
	 * @return knot vector
	 */
	public double[] getKnots() {
		return knots;
	}

	/**
	 * Set new knot vector for the NURBS curve.
	 * 
	 * @param newKnots
	 *            new nonperiodic knot vector
	 */
	public void setKnots(double[] newKnots) {
		// need checking here!
		int n = controlPoints.length - 1;

		if (!Algorithms.IsKnotsValid(n, degree, newKnots))
			throw new GeomException(
					"setKnots() failed! Not a valid knot vector!");
		knots = newKnots;
	}

	public String toString() {
		String s;
		s = "NURBS Curve\n";
		s = s + "  Degree: " + Integer.toString(degree) + "\n";
		s = s + "  Number of ControlPoints: "
				+ Integer.toString(controlPoints.length) + "\n";
		for (int i = 0; i < controlPoints.length; i++) {
			s = s + "    " + controlPoints[i].toString() + "\n";
		}
		return s;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see fq.geom.ParametricCurve#getParamExtents()
	 */
	public double[] getParamExtents() {
		double[] r = { this.minKnot(), this.maxKnot() };
		return r;
	}

	public int EntityType() {
		return 126;
	}
}
