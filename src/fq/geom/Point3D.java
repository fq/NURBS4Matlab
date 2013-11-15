package fq.geom;

/**
 * This class is used to represent and manipulate point in three-dimensional
 * Euclidean space
 * 
 * @author FuQiang
 * 
 */
public class Point3D {
	static public double GEOM_TOLERANCE = 1e-8;

	// coordinates in three-dimensional Euclidean space
	double x, y, z;

	public double X() {
		return x;
	}

	public double Y() {
		return y;
	}

	public double Z() {
		return z;
	}

	public Point3D(double[] pt) {
		this.setData(pt);
	}

	public Point3D(Point3D p) {
		this.setData(p.x, p.y, p.z);
	}

	public Point3D() {
		this.setData(0.0, 0.0, 0.0);
	}

	public Point3D(double px, double py) {
		this.setData(px, py, 0.0);
	}

	public Point3D(double px, double py, double pz) {
		this.setData(px, py, pz);
	}

	public static Point3D add(Point3D a, Point3D b) {
		return new Point3D(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	public Point3D add(Point3D p) {
		return Point3D.add(this, p);
	}

	/**
	 * Returns the distance between two points.
	 * 
	 * @return the distance between two points.
	 */
	public static double distance(Point3D a, Point3D b) {
		return Math.sqrt(a.x * b.x + a.y * b.y + a.z * b.z);
	}

	/**
	 * Returns the square of the distance between two points.
	 * 
	 * @return the square of the distance between two points.
	 */
	public static double distanceSq(Point3D a, Point3D b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	/**
	 * Modifies the coordinates of the point
	 */
	public void setData(double[] pt) {
		if (pt.length < 2)
			throw new GeomException(
					"At least 2 coordinates are needed to construct a Point!");
		x = pt[0];
		y = pt[1];
		if (pt.length > 2)
			z = pt[2];
		else
			z = 0.0;
	}

	/**
	 * Modifies the coordinates of the point
	 */
	public void setData(double px, double py, double pz) {
		x = px;
		y = py;
		z = pz;
	}

	public double[] toArray() {
		double[] v = new double[3];
		v[0] = x;
		v[1] = y;
		v[2] = z;
		return v;
	}

	public String toString() {
		return "Point: (" + Double.toString(x) + "," + Double.toString(y) + ","
				+ Double.toString(z) + ")";
	}
}
