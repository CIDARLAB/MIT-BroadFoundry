import java.util.ArrayList;

public class ScratchWork{
	public static void main(String[] args){
		ConstantProperties.setOpNames();
		String aName = "(((IN2.IN3).(IN3.IN4))+(((IN3.IN4).IN1).IN2))";
		String bName = "(~((((~IN2).IN3).((~IN2).IN4)).(IN3.IN4)))";
		String cName = "(~(((IN2.IN4).IN1).(IN2.IN3)))";
		String dName = "(~(((((IN3.IN4).IN3).((IN3.IN4).IN4)).(~IN2)).(((IN3.IN4).IN3).IN2)))";
		String eName = "(((~IN2).IN3).IN4)";
		String fName = "(~((((~IN3).IN2).IN4).((IN1.IN2).IN3)))";
		String gName = "(~((((~IN3).(IN1.IN4)).(IN1.IN2)).((~IN3).IN2)))";
		Circuit a = new Circuit(aName,"0");
		Circuit b = new Circuit(bName,"0");
		Circuit c = new Circuit(cName,"0");
		Circuit d = new Circuit(dName,"0");
		Circuit e = new Circuit(eName,"0");
		Circuit f = new Circuit(fName,"0");
		Circuit g = new Circuit(gName,"0");
		System.out.println("a: "+a.getNetlist());
		System.out.println("b: "+b.getNetlist());
		System.out.println("c: "+c.getNetlist());
		System.out.println("d: "+d.getNetlist());
		System.out.println("e: "+e.getNetlist());
		System.out.println("f: "+f.getNetlist());
		System.out.println("g: "+g.getNetlist());
	}
}