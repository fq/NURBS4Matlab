package fq.geom;

/**
 * Abstract parametric surface class<br>
 * This class is used as a basis for all parametric surfaces.
 * 
 * @author FuQiang
 * 
 */
public abstract class ParametricSurface {	
	/**
	 * Compuate a point on the surface at parametric value (u,v)
	 * 
	 * @param u -
	 *            the parametric value
	 * @param v -
	 *            the parametric value
	 */
	public abstract Point3D getPointAt(double u, double v);

	/**
	 * Compuate derivates of the surface at parametric value (u,v)
	 * 
	 * @param u -
	 *            the u parametric value
	 * @param v -
	 *            the v parametric value
	 * @param d -
	 *            the derivative is computed up to and including to this value
	 * 
	 * @return An 2-dimensional array of Vector3D containing the derivatives<br>
	 *         SKL(k,l), while 0 <= k+l <= d :<br>
	 *         &nbsp k - the derivative with respect to <i>u k</i> times<br>
	 *         &nbsp l - the derivative with respect to <i>v l</i> times
	 */
	public abstract Vector3D[][] getDerivsAt(double u, double v, int d);

	/**
	 * Returns the parametric extents in U-direction of the surface.
	 * 
	 * @return params[0] - The parameter associated with the start point of the
	 *         U-direction<br>
	 *         params[1] - The parameter associated with the end point of the
	 *         U-direction
	 */
	public abstract double[] getParamExtentsU();

	/**
	 * Returns the parametric extents in V-direction of the surface.
	 * 
	 * @return params[0] - The parameter associated with the start point of the
	 *         V-direction<br>
	 *         params[1] - The parameter associated with the end point of the
	 *         V-direction
	 */
	public abstract double[] getParamExtentsV();

	/**
	 * Returns the parametric extents in U-direction of the surface.
	 * (This method is designed for using in MATLAB)
	 * 
	 * @return params[0] - The parameter associated with the start point of the
	 *         U-direction<br>
	 *         params[1] - The parameter associated with the end point of the
	 *         U-direction
	 */
	public double[] ParamExtentsU()
	{
		return this.getParamExtentsU();
	}

	/**
	 * Returns the parametric extents in V-direction of the surface.
	 * (This method is designed for using in MATLAB)
	 * 
	 * @return params[0] - The parameter associated with the start point of the
	 *         V-direction<br>
	 *         params[1] - The parameter associated with the end point of the
	 *         V-direction
	 */
	public double[] ParamExtentsV()
	{
		return this.getParamExtentsV();
	}


	/**
	 * Returns the unit Normal Vector at the parameter specified
	 * 
	 * @param u -
	 *            the parametric value
	 * @param v -
	 *            the parametric value
	 */
	public Vector3D getNormalAt(double u, double v) {
		Vector3D[][] ders = this.getDerivsAt(u, v, 1);
		Vector3D du = ders[1][0];
		Vector3D dv = ders[0][1];

		Vector3D r = Vector3D.crossProduct(du, dv);
		if (r.Magnitude() < Vector3D.GEOM_TOLERANCE) {
			// cannot calcuate normal using u and v direction
			// System.out.println(r.Magnitude());
			return new Vector3D(0.0, 0.0, 0.0);
		} else {
			r.Normalize();
			return r;
		}
	}

	/**
	 * Compuate a point on the surface at a parametric value (u,v)<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parametric value
	 * @param v -
	 *            the parametric value
	 * @return the coordinates of the point, p[3]
	 */
	public double[] PointAt(double u, double v) {
		double[] paramU = this.getParamExtentsU();
		if (u < paramU[0])
			u = paramU[0];
		if (u > paramU[1])
			u = paramU[1];
		double[] paramV = this.getParamExtentsV();
		if (v < paramV[0])
			v = paramV[0];
		if (v > paramV[1])
			v = paramV[1];
		Point3D p = this.getPointAt(u, v);
		return p.toArray();
	}

	/**
	 * Generate a surface mesh from the input parameters in U & V directions<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            array of parameter values, U direction
	 * @param v -
	 *            array of parameter values, V direction
	 * @return array of points on the surface, p(NU,NV,3)
	 */
	public double[][][] PointAt(double[] u, double[] v) {
		int nu = u.length;
		int nv = v.length;

		double[][][] p = new double[nu][nv][];
		for (int i = 0; i < nu; i++)
			for (int j = 0; j < nv; j++) {
				p[i][j] = this.PointAt(u[i], v[j]);
			}
		return p;
	}

	/**
	 * Generate N points from the input parameters <b>uv(2,N)</b><br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param uv -
	 *            uv(0,:) is parameter values in U direction; uv(1,:) is
	 *            parameter values in V direction
	 * @return array of points on the curve, p(NU,NV,3)
	 */
	public double[][] PointAt(double[][] uv) {
		int n = uv[0].length;

		double[][] p = new double[n][];
		for (int i = 0; i < n; i++) {
			p[i] = this.PointAt(uv[0][i], uv[1][i]);
		}
		return p;
	}

	/**
	 * Compuate derivates of the surface at parametric value (u,v)<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the u parametric value
	 * @param v -
	 *            the v parametric value
	 * @param d -
	 *            the derivative is computed up to and including to this value
	 * 
	 * @return An 3-dimensional array SKL containing the derivatives<br>
	 *         SKL(k,l,3), while 0 <= k+l <= d :<br>
	 *         &nbsp k - the derivative with respect to <i>u k</i> times<br>
	 *         &nbsp l - the derivative with respect to <i>v l</i> times<br>
	 *         &nbsp 3 - derivative values, x y z<br>
	 */
	public double[][][] DerivsAt(double u, double v, int d) {
		double[] paramU = this.getParamExtentsU();
		if (u < paramU[0])
			u = paramU[0];
		if (u > paramU[1])
			u = paramU[1];
		double[] paramV = this.getParamExtentsV();
		if (v < paramV[0])
			v = paramV[0];
		if (v > paramV[1])
			v = paramV[1];

		Vector3D[][] ders = this.getDerivsAt(u, v, d);
		double[][][] skl = new double[d + 1][d + 1][];
		for (int i = 0; i <= d; i++)
			for (int j = 0; j <= d; j++) {
				skl[i][j] = ders[i][j].toArray();
			}
		return skl;
	}

	/**
	 * Returns the unit Normal Vector at the parameter specified<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            the parametric value
	 * @param v -
	 *            the parametric value
	 */
	public double[] NormalAt(double u, double v) {
		double[] paramU = this.getParamExtentsU();
		if (u < paramU[0])
			u = paramU[0];
		if (u > paramU[1])
			u = paramU[1];
		double[] paramV = this.getParamExtentsV();
		if (v < paramV[0])
			v = paramV[0];
		if (v > paramV[1])
			v = paramV[1];

		Vector3D r = this.getNormalAt(u, v);
		return r.toArray();
	}

	/**
	 * Generate a set of Normal vectors on the surface<br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param u -
	 *            array of parameter values, U direction
	 * @param v -
	 *            array of parameter values, V direction
	 * @return array of points on the curve, p(NU,NV,3)
	 */
	public double[][][] NormalAt(double[] u, double[] v) {
		int nu = u.length;
		int nv = v.length;

		double[][][] p = new double[nu][nv][];
		for (int i = 0; i < nu; i++)
			for (int j = 0; j < nv; j++) {
				p[i][j] = getNormalAt(u[i], v[j]).toArray();
			}
		return p;
	}
	
	/**
	 * Generate N normal vector from the input parameters <b>uv(2,N)</b><br>
	 * (This method is designed for using in MATLAB)
	 * 
	 * @param uv -
	 *            uv(0,:) is parameter values in U direction; uv(1,:) is
	 *            parameter values in V direction
	 * @return array of normal vector on the curve, p(NU,NV,3)
	 */
	public double[][] NormalAt(double[][] uv) {
		int n = uv[0].length;

		double[][] p = new double[n][];
		for (int i = 0; i < n; i++) {
			p[i] = this.NormalAt(uv[0][i], uv[1][i]);
		}
		return p;
	}

	public int EntityType(){
		return 128;
	}
}
