package fq.geom;

/**
 * Abstract parametric curve class<br>
 * This class is used as a basis for all parametric curves.
 * 
 * @author FuQiang
 * 
 */
public abstract class ParametricCurve {
	/**
	 * Compute dth derivates of the curve at a parametric value
	 * 
	 * @param u -
	 *            the parametric value
	 * @param d -
	 *            the derivative is computed up to and including to this value
	 * @return An Vector3D array of derivatives, Ders[k], k from 0 to d<br>
	 *         Ders[0] -> the point<br>
	 *         Ders[1] -> the 1th derivative<br>
	 *         ......<br>
	 *         Ders[k] -> the kth derivative.
	 */
	public abstract Vector3D[] getDerivsAt(double u, int d);

	/**
	 * Compute a point on the curve at a parametric value
	 * 
	 * @param u -
	 *            the parametric value
	 * @return a point on the curve, Point3D
	 */
	public abstract Point3D getPointAt(double u);

	/**
	 * Returns the parametric extents of the curve. This is the parametric
	 * equivalent of the end-points.
	 * 
	 * @return params[0] - The parameter associated with the start point of the
	 *         curve<br>
	 *         params[1] - The parameter associated with the end point of the
	 *         curve
	 */
	public abstract double[] getParamExtents();

	/**
	 * Returns the parametric extents of the curve. This is the parametric
	 * equivalent of the end-points.
	 * (This method is designed for using in MATLAB)
	 * 
	 * @return params[0] - The parameter associated with the start point of the
	 *         curve<br>
	 *         params[1] - The parameter associated with the end point of the
	 *         curve
	 */
	public double[] ParamExtents()
	{
		return this.getParamExtents();
	}

	/**
	 * Returns the absolute curvature at the parameter specified
	 * 
	 * @param u -
	 *            the parametric value
	 */
	public double getCurvatureAt(double u) {
		Vector3D[] ders = this.getDerivsAt(u, 2);
		Vector3D p1 = ders[1];
		Vector3D p2 = ders[2];

		Vector3D p = Vector3D.crossProduct(p1, p2);
		double t = p1.Magnitude();
		return p.Magnitude() / (t * t * t);
	}

	/**
	 * Returns the unit-vector tangent at the parameter specified
	 * 
	 * @param u -
	 *            the parametric value
	 */
	public Vector3D getTangentAt(double u) {
		Vector3D[] ders = this.getDerivsAt(u, 1);
		Vector3D v = ders[1];
		v.Normalize();
		return v;
	}

	/**
	 * Returns the unit principal normal vector at the parameter specified
	 * 
	 * @param u -
	 *            the parametric value
	 */
	public Vector3D getNormalAt(double u) {
		Vector3D[] ders = this.getDerivsAt(u, 2);
		Vector3D p1 = ders[1];
		Vector3D p2 = ders[2];

		Vector3D t = ders[1];
		t.Normalize();

		Vector3D b = Vector3D.crossProduct(p1, p2);
		b.Normalize();

		Vector3D n = Vector3D.crossProduct(b, t);
		return n;
	}

	/**
	 * Compuate dth derivates of the curve at a parametric value.<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parametric value
	 * @param d -
	 *            the derivative is computed up to and including to this value
	 * @return An array of the kth derivative, Ders[k][3], k from 0 to d<br>
	 *         Ders[0][] -> the point<br>
	 *         Ders[1][] -> the 1th derivative<br>
	 *         ......<br>
	 *         Ders[k][] -> the kth derivative.
	 */
	public double[][] DerivsAt(double u, int d) {
		double[] param = this.getParamExtents();
		if (u < param[0])
			u = param[0];
		if (u > param[1])
			u = param[1];
		Vector3D[] ders = this.getDerivsAt(u, d);
		double[][] r = new double[d + 1][];
		for (int i = 0; i <= d; i++) {
			r[i] = ders[i].toArray();
		}
		return r;
	}

	/**
	 * Compuate a point on the curve at a parametric value<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parameter value
	 * @return the coordinates of the point, p[3]
	 */
	public double[] PointAt(double u) {
		double[] param = this.getParamExtents();
		if (u < param[0])
			u = param[0];
		if (u > param[1])
			u = param[1];
		return this.getPointAt(u).toArray();
	}

	/**
	 * Generate a set of points on the curve<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            array of parameter values
	 * @return array of points on the curve, p(n,3)
	 */
	public double[][] PointAt(double[] u) {
		int n = u.length;

		double[][] p = new double[n][];
		for (int i = 0; i < n; i++) {
			p[i] = this.PointAt(u[i]);
		}
		return p;
	}

	/**
	 * Returns the unit tangent vector at the parameter specified<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parametric value
	 */
	public double[] TangentAt(double u) {
		double[] param = this.getParamExtents();
		if (u < param[0])
			u = param[0];
		if (u > param[1])
			u = param[1];
		return this.getTangentAt(u).toArray();
	}

	/**
	 * Generate a set of unit tangent vectors on the curve<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            array of parameter values
	 * @return array of unit tangent vector on the curve, t(n,3)
	 */
	public double[][] TangentAt(double[] u) {
		int n = u.length;

		double[][] t = new double[n][];
		for (int i = 0; i < n; i++) {
			t[i] = this.TangentAt(u[i]);
		}
		return t;
	}

	/**
	 * Returns the curvature at the parameter specified<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parametric value
	 */
	public double CurvatureAt(double u) {
		double[] param = this.getParamExtents();
		if (u < param[0])
			u = param[0];
		if (u > param[1])
			u = param[1];
		return this.getCurvatureAt(u);
	}

	/**
	 * Returns a set of curvature at the parameter specified<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            array of parameter values
	 */
	public double[] CurvatureAt(double[] u) {
		int n = u.length;

		double[] ck = new double[n];
		for (int i = 0; i < n; i++) {
			ck[i] = this.CurvatureAt(u[i]);
		}
		return ck;
	}

	/**
	 * Returns the unit proncipal normal vector at the parameter specified<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parametric value
	 */
	public double[] NormalAt(double u) {
		double[] param = this.getParamExtents();
		if (u < param[0])
			u = param[0];
		if (u > param[1])
			u = param[1];
		return this.getNormalAt(u).toArray();
	}

	/**
	 * Returns a set of unit principal normal vectors<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parametric value
	 * @return array of unit proncipal normal vectors on the curve, N(n,3)
	 */
	public double[][] NormalAt(double[] u) {
		int n = u.length;

		double[][] t = new double[n][];
		for (int i = 0; i < n; i++) {
			t[i] = this.NormalAt(u[i]);
		}
		return t;
	}
}
