package fq.geom.nurbs;

import fq.geom.*;

/**
 * NURBS surface class<br>
 * This class is used to represent and manipulate NURBS surface.
 * 
 * @author FuQiang
 * 
 */
public class NurbsSurface extends ParametricSurface {

	double[] knotsU; // knot vector - U direction

	double[] knotsV; // knot vector - V direction

	Point3DH[][] controlPoints; // control points

	int degreeU, degreeV; // degrees

	public NurbsSurface(double[][][] cp) {
		init(cp, 3, 3);
	}

	/**
	 * Initialize the NURRBS surface
	 * 
	 * @param cp -
	 *            3-D array of (NU x NV) control points, cp(NU,NV,4)<br>
	 *            cp(i,j,1) -> x coordinates of the (i,j) control points<br>
	 *            cp(i,j,2) -> y coordinates of the (i,j) control points<br>
	 *            cp(i,j,3) -> z coordinates of the (i,j) control
	 *            points(optional, default is 0.0)<br>
	 *            cp(i,j,4) -> weights of the (i,j) control points(optional,
	 *            default is 1.0)<br>
	 *            where, 0 <= i <= NU, 0 <= j <= NV;
	 * @param degU -
	 *            Degree of the NURBS surface at U-direction
	 * @param degV -
	 *            Degree of the NURBS surface at V-direction
	 */
	public NurbsSurface(double[][][] cp, int degU, int degV) {
		init(cp, degU, degV);
	}

	/**
	 * Initialize the NURBS surface
	 */
	void init(double[][][] cp, int degU, int degV) {
		int nu = cp.length - 1;
		int nv = cp[0].length - 1;
		if (degU < 1 || degV < 1)
			throw new GeomException("degree must >= 1!");
		if (nu < 1 || nv < 1)
			throw new GeomException("at least two control points!");
		if (degU >= nu)
			degU = nu;
		if (degV >= nv)
			degV = nv;

		controlPoints = new Point3DH[nu + 1][nv + 1];
		for (int i = 0; i <= nu; i++)
			for (int j = 0; j <= nv; j++) {
				controlPoints[i][j] = new Point3DH(cp[i][j]);
			}
		degreeU = degU;
		degreeV = degV;
		knotsU = Algorithms.GenKnots(nu, degreeU);
		knotsV = Algorithms.GenKnots(nv, degreeV);
	}

	/**
	 * Determine the U direction knot span index at u
	 * 
	 * @param u -
	 *            the parametric value, U-direction
	 */
	public int findSpanU(double u) {
		return Algorithms.FindSpan(u, degreeU, knotsU);
	}

	/**
	 * Determine the V direction knot span index at v
	 * 
	 * @param v -
	 *            the parametric value, V-direction
	 */
	public int findSpanV(double v) {
		return Algorithms.FindSpan(v, degreeV, knotsV);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see fq.geom.ParametricSurface#getPointAt(double, double)
	 */
	public Point3D getPointAt(double u, double v) {
		if (u < minKnotU() || u > maxKnotU())
			throw new GeomException("Input parameter out of range!");
		if (v < minKnotV() || v > maxKnotV())
			throw new GeomException("Input parameter out of range!");

		int uspan = findSpanU(u);
		int vspan = findSpanV(v);
		double[] Nu = Algorithms.BasisFuns(u, uspan, degreeU, knotsU);
		double[] Nv = Algorithms.BasisFuns(v, vspan, degreeV, knotsV);
		double[][] temp = new double[degreeV + 1][];

		int l;
		for (l = 0; l <= degreeV; l++) {
			temp[l] = new double[4];
			for (int k = 0; k <= degreeU; k++) {
				temp[l][0] += Nu[k]
						* controlPoints[uspan - degreeU + k][vspan - degreeV
								+ l].XW();
				temp[l][1] += Nu[k]
						* controlPoints[uspan - degreeU + k][vspan - degreeV
								+ l].YW();
				temp[l][2] += Nu[k]
						* controlPoints[uspan - degreeU + k][vspan - degreeV
								+ l].ZW();
				temp[l][3] += Nu[k]
						* controlPoints[uspan - degreeU + k][vspan - degreeV
								+ l].W();
			}
		}
		double xw = 0.0, yw = 0.0, zw = 0.0, ww = 0.0;
		for (l = 0; l <= degreeV; l++) {
			xw += Nv[l] * temp[l][0];
			yw += Nv[l] * temp[l][1];
			zw += Nv[l] * temp[l][2];
			ww += Nv[l] * temp[l][3];
		}
		return new Point3D(xw / ww, yw / ww, zw / ww);
	}

	/**
	 * Compuate derivates on the NURBS surface at homogeneous coordinates
	 */
	private double[][][] getDerivsAtH(double u, double v, int d) {
		if (u < minKnotU() || u > maxKnotU())
			throw new GeomException("Input parameter out of range!");
		if (v < minKnotV() || v > maxKnotV())
			throw new GeomException("Input parameter out of range!");

		int k, l, du, dv;
		double[][][] skl = new double[d + 1][d + 1][4];

		du = Math.min(d, degreeU);
		dv = Math.min(d, degreeV);

		int uspan = findSpanU(u);
		int vspan = findSpanV(v);
		double[][] Nu = Algorithms.DersBasisFuns(u, du, uspan, degreeU, knotsU);
		double[][] Nv = Algorithms.DersBasisFuns(v, dv, vspan, degreeV, knotsV);

		double[][] temp = new double[degreeV + 1][];
		int dd, r, s;
		for (k = 0; k <= du; k++) {
			for (s = 0; s <= degreeV; s++) {
				temp[s] = new double[4];
				for (r = 0; r <= degreeU; r++) {
					temp[s][0] += Nu[k][r]
							* controlPoints[uspan - degreeU + r][vspan
									- degreeV + s].XW();
					temp[s][1] += Nu[k][r]
							* controlPoints[uspan - degreeU + r][vspan
									- degreeV + s].YW();
					temp[s][2] += Nu[k][r]
							* controlPoints[uspan - degreeU + r][vspan
									- degreeV + s].ZW();
					temp[s][3] += Nu[k][r]
							* controlPoints[uspan - degreeU + r][vspan
									- degreeV + s].W();
				}
			}
			dd = Math.min(d - k, dv);
			for (l = 0; l <= dd; l++) {
				// skl[k][l] = new double[4];
				skl[k][l][0] = 0.0;
				skl[k][l][1] = 0.0;
				skl[k][l][2] = 0.0;
				skl[k][l][3] = 0.0;
				for (s = 0; s <= degreeV; s++) {
					skl[k][l][0] += Nv[l][s] * temp[s][0];
					skl[k][l][1] += Nv[l][s] * temp[s][1];
					skl[k][l][2] += Nv[l][s] * temp[s][2];
					skl[k][l][3] += Nv[l][s] * temp[s][3];
				}
			}
		}
		return skl;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see fq.geom.ParametricSurface#getDerivsAt(double, double, int)
	 */
	public Vector3D[][] getDerivsAt(double u, double v, int d) {
		if (u < minKnotU() || u > maxKnotU())
			throw new GeomException("Input parameter out of range!");
		if (v < minKnotV() || v > maxKnotV())
			throw new GeomException("Input parameter out of range!");

		int k, l;
		// double[][][] skl = new double[d + 1][d + 1][3];
		Vector3D[][] skl = new Vector3D[d + 1][d + 1];
		int i, j;
		for (i = 0; i <= d; i++)
			for (j = 0; j <= d; j++) {
				skl[i][j] = new Vector3D();
			}

		double[][][] ders = getDerivsAtH(u, v, d);
		double[][] Bin = Algorithms.BinomialCoef(d, d);
		double[] pv = new double[3];
		double[] pv2;

		for (k = 0; k <= d; k++) {
			for (l = 0; l <= d - k; l++) {
				pv[0] = ders[k][l][0];
				pv[1] = ders[k][l][1];
				pv[2] = ders[k][l][2];
				for (j = 1; j <= l; j++) {
					pv[0] -= Bin[l][j] * ders[0][j][3] * skl[k][l - j].X();
					pv[1] -= Bin[l][j] * ders[0][j][3] * skl[k][l - j].Y();
					pv[2] -= Bin[l][j] * ders[0][j][3] * skl[k][l - j].Z();
				}
				for (i = 1; i <= k; i++) {
					pv[0] -= Bin[k][i] * ders[i][0][3] * skl[k - i][l].X();
					pv[1] -= Bin[k][i] * ders[i][0][3] * skl[k - i][l].Y();
					pv[2] -= Bin[k][i] * ders[i][0][3] * skl[k - i][l].Z();
					pv2 = new double[3];
					for (j = 1; j <= l; j++) {
						pv2[0] += Bin[l][j] * ders[i][j][3]
								* skl[k - i][l - j].X();
						pv2[1] += Bin[l][j] * ders[i][j][3]
								* skl[k - i][l - j].Y();
						pv2[2] += Bin[l][j] * ders[i][j][3]
								* skl[k - i][l - j].Z();
					}
					pv[0] -= Bin[k][i] * pv2[0];
					pv[1] -= Bin[k][i] * pv2[1];
					pv[2] -= Bin[k][i] * pv2[2];
				}
				skl[k][l].setData(pv[0] / ders[0][0][3], pv[1] / ders[0][0][3],
						pv[2] / ders[0][0][3]);
			}
		}
		return skl;
	}

	public double minKnotU() {
		return knotsU[0];
	}

	public double maxKnotU() {
		return knotsU[knotsU.length - 1];
	}

	public double minKnotV() {
		return knotsV[0];
	}

	public double maxKnotV() {
		return knotsV[knotsV.length - 1];
	}

	/**
	 * Returns the control points of the NURBS surface
	 * 
	 * @return An 3-dimensional array CP of control points<br>
	 *         CP(U,V,4):<br>
	 *         &nbspU - Number of control points in U direction.<br>
	 *         &nbspV - Number of control points in V direction.<br>
	 *         &nbsp4 - Coordinates of control point(U,V), including weights.<br>
	 * 
	 */
	public double[][][] getControlPoints() {
		double[][][] p = new double[controlPoints.length][controlPoints[0].length][];
		for (int i = 0; i < controlPoints.length; i++)
			for (int j = 0; j < controlPoints[0].length; j++) {
				p[i][j] = controlPoints[i][j].toArray();
			}
		return p;
	}

	/**
	 * Returns the degree in U-direction
	 */
	public int getDegreeU() {
		return degreeU;
	}

	/**
	 * Returns the degree in V-direction
	 */
	public int getDegreeV() {
		return degreeV;
	}

	/**
	 * Returns the knot vector in U-direction
	 * 
	 * @return knot vector in U-direction
	 */
	public double[] getKnotsU() {
		return knotsU;
	}

	/**
	 * Set new knot vector in U-direction
	 * 
	 * @param newKnots
	 *            new nonperiodic knot vector
	 */
	public void setKnotsU(double[] newKnots) {
		int nu = controlPoints.length - 1;
		if (!Algorithms.IsKnotsValid(nu, degreeU, newKnots))
			throw new GeomException(
					"setKnotsU() failed! Not a valid knot vector!");
		knotsU = newKnots;
	}

	/**
	 * Set new knot vector in V-direction
	 * 
	 * @param newKnots
	 *            new nonperiodic knot vector
	 */
	public void setKnotsV(double[] newKnots) {
		int nv = controlPoints[0].length - 1;
		if (!Algorithms.IsKnotsValid(nv, degreeV, newKnots))
			throw new GeomException(
					"setKnotsV() failed! Not a valid knot vector!");
		knotsV = newKnots;
	}

	/**
	 * Returns the knot vector in V-direction
	 * 
	 * @return knot vector in V-direction
	 */
	public double[] getKnotsV() {
		return knotsV;
	}

	public String toString() {
		String s;
		int nu = controlPoints.length;
		int nv = controlPoints[0].length;
		s = "NURBS Surface\n";
		s = s + "  Degree: (" + Integer.toString(degreeU) + " x "
				+ Integer.toString(degreeV) + ")\n";
		s = s + "  Number of ControlPoints: (" + Integer.toString(nu) + " x "
				+ Integer.toString(nv) + ")\n";
		return s;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see fq.geom.ParametricSurface#GetParamExtentsU()
	 */
	public double[] getParamExtentsU() {
		double[] r = { this.minKnotU(), this.maxKnotU() };
		return r;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see fq.geom.ParametricSurface#GetParamExtentsV()
	 */
	public double[] getParamExtentsV() {
		double[] r = { this.minKnotV(), this.maxKnotV() };
		return r;
	}

	public double[] EvalSP1(double u, double v){
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
		
		double[] duv = new double[2];
		Vector3D[][] ders = this.getDerivsAt(u, v, 1);
		Vector3D du = ders[1][0];
		Vector3D dv = ders[0][1];

		Vector3D n = Vector3D.crossProduct(du, dv);
		n.Normalize();
		double tmp = Vector3D.dotProduct(n,new Vector3D(0,0,1));
		double spx,spy;
		spx = - n.X() * tmp;
		spy = - n.Y() * tmp;
		
		duv[0] = (dv.Y()*spx - dv.X()*spy);
		duv[1] = - (du.Y()*spx - du.X()*spy);
		
		return duv;
	}
	
	public double[] EvalSP2(double u,double v){
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
		
		double[] duv = new double[2];

		int uspan = findSpanU(u);
		int vspan = findSpanV(v);
		double[][] U = Algorithms.DersBasisFuns(u, 1, uspan, degreeU, knotsU);
		double[][] V = Algorithms.DersBasisFuns(v, 1, vspan, degreeV, knotsV);
		
		int s,r;
		double xw = 0.0,yw = 0.0,zw =0.0,ww=0.0;
		double xw_u = 0.0,yw_u = 0.0,zw_u =0.0,ww_u=0.0;
		double xw_v = 0.0,yw_v = 0.0,zw_v =0.0,ww_v=0.0;
		double[] tmp;
		for (s = 0; s<=degreeV; s++){
			tmp = new double[8];
			for(r = 0;r <= degreeU; r++) {
				tmp[0] += U[0][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].XW();
				tmp[1] += U[0][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].YW();
				tmp[2] += U[0][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].ZW();
				tmp[3] += U[0][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].W();
				
				tmp[4] += U[1][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].XW();
				tmp[5] += U[1][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].YW();
				tmp[6] += U[1][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].ZW();
				tmp[7] += U[1][r] * controlPoints[uspan - degreeU + r][vspan - degreeV + s].W();
			}
			xw += tmp[0] * V[0][s];
			yw += tmp[1] * V[0][s];
			zw += tmp[2] * V[0][s];
			ww += tmp[3] * V[0][s];
			
			xw_u += tmp[4] * V[0][s];
			yw_u += tmp[5] * V[0][s];
			zw_u += tmp[6] * V[0][s];
			ww_u += tmp[7] * V[0][s];

			xw_v += tmp[0] * V[1][s];
			yw_v += tmp[1] * V[1][s];
			zw_v += tmp[2] * V[1][s];
			ww_v += tmp[3] * V[1][s];
		}
		double[] M1 = {
				xw * ww_v - xw_v * ww,
				yw * ww_v - yw_v * ww};
		double a1,a2,a3;
		a1 = xw * M1[0] + yw * M1[1];
		a2 = xw_u * M1[0] + yw_u * M1[1];
		a3 = xw_v * M1[0] + yw_v * M1[1];
		
		double[] M2 = {
				- xw * ww_u + xw_u * ww,
				- yw * ww_u + yw_u * ww} ;
		double b1,b2,b3;
		b1 = xw * M2[0] + yw * M2[1];
		b2 = xw_u * M2[0] + yw_u * M2[1];
		b3 = xw_v * M2[0] + yw_v * M2[1];
		
		duv[0] = a1*zw_u*ww_v + zw*ww_u*a3 + ww*a2*zw_v
				-ww*zw_u*a3 - a1*ww_u*zw_v - zw*a2*ww_v;
		duv[1] = b1*zw_u*ww_v + zw*ww_u*b3 + ww*b2*zw_v
				-ww*zw_u*b3 - b1*ww_u*zw_v - zw*b2*ww_v;
		
		return duv;
	}
}
