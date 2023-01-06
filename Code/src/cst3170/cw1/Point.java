package cst3170.cw1;

public class Point {

	double x;
	double y;
	
//zero argument constructor	
	public Point() {
		this.x=0;
		this.y=0;
	}
//parameterised Constructor	
	public Point(double a, double b)
	{
		this.x=a;
		this.y=b;
	}
//getters	
	public double getX() {
		return this.x;
	}
	public double getY() {
		return this.y;
	}
//setters	
	public void setX(double a) {
		this.x = a;
	}
	public void setY(double a) {
		this.y = a;
	}
	

}
