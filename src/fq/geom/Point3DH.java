package fq.geom;

/**
 * 	This class is used to represent and manipulate points in 
 * 	homogenous coordinates.
 * 
 * @author FuQiang
 *
 */
public class Point3DH extends Point3D {
	// homogenous coordiante
	double w;
	
	public double XW() {
		return x * w;
	}

	public double YW() {
		return y * w;
	}

	public double ZW() {
		return z * w;
	}

	public double W() {
		return w;
	}
	
	public Point3DH() {
		this(0.0,0.0,0.0,1.0);
	}
	
	public Point3DH(double[] pt) {
		super(pt);
		if (pt.length > 3) {
//			if (pt[3] == 0)
	//			throw new GeomException("Weights cannot equal to 0!");
			w = pt[3];
		} else
			w = 1.0;
	}
	
	public Point3DH(double px, double py, double pz, double pw) {
		super(px,py,pz);
		//if (pw == 0)
			//throw new GeomException("Weights cannot equal to 0!");
		w = pw;
	}
	
	public double[] toArray() {
		double[] v = new double[4];
		v[0] = x;
		v[1] = y;
		v[2] = z;
		v[3] = w;
		return v;
	}
	
	public String toString() {
		return "PointH: ("  + Double.toString(x) + "," + Double.toString(y)
		+ "," + Double.toString(z) + "," + Double.toString(w) + ")";
	}
}
