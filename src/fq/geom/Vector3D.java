package fq.geom;

public class Vector3D extends Point3D {

	public Vector3D(){
		super();
	}

	public Vector3D(double[] pt){
		super(pt);
	}
	
	public Vector3D(double px,double py){
		super(px,py);
	}

	public Vector3D(double px,double py,double pz){
		super(px,py,pz);
	}
	
	public Vector3D(Vector3D p) {
		this.setData(p.x, p.y, p.z);
	}
	
	public Vector3D(Point3D p) {
		this.setData(p.x, p.y, p.z);
	}

	public static Vector3D add(Vector3D a, Vector3D b){
		return new Vector3D(a.x+b.x, a.y+b.y, a.z+b.z);
	}

	public Vector3D add(Vector3D p) {
		return Vector3D.add(this,p);
	}

	public static Vector3D crossProduct(Vector3D a, Vector3D b) {
		Vector3D r = new Vector3D();
		r.x = a.y * b.z - a.z * b.y;
		r.y = a.z * b.x - a.x * b.z;
		r.z = a.x * b.y - a.y * b.x;
		return r;
	}

	public static double dotProduct(Vector3D a, Vector3D b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	
	public Vector3D crossProduct(Vector3D v) {
		return Vector3D.crossProduct(this, v);
	}

	public double dotProduct(Vector3D v) {
		return Vector3D.dotProduct(this, v);
	}
	
	public void Normalize() {
		double l = Magnitude();
		if (l != 0.0) {
			x = x / l;
			y = y / l;
			z = z / l;
		}
	}
	
	public double Magnitude() {
		return Math.sqrt(x * x + y * y + z * z);
	}

	public String toString() {
		return "Vector3D: (" + Double.toString(x) + "," + Double.toString(y) + ","
				+ Double.toString(z) + ")";
	}
}
