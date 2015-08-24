import org.python.util.PythonInterpreter;
import org.python.core.*;
//Download Jython 2.7.0 - Standalone Jar from:
//http://www.jython.org/downloads.html
//Alternatively you can get the jar file from the project folder.

public class Practice{
	public static void main(String[] args){
		System.out.println("start");
		double startTime1 = System.currentTimeMillis();
		
		PythonInterpreter python = new PythonInterpreter();
		double startTime2 = System.currentTimeMillis();
		System.out.println("start1");
		for (int i=0;i<1;i++){
			/*
			//System.out.println("running");
			//You can use absolute path or just put the name of the file if the file is in the project folder
			python.execfile("C:/Users/Arinze/OneDrive/UROP_Summer_2015/RunPythonInJava/practiceCode.py");
			python.execfile("practiceCode.py");
			python.exec("he = func()");
			PyObject he1 = python.get("he");
			//System.out.println(he1);
			String h = he1.asString();
			//System.out.println(h);
			
			int number1 = 10;
			int number2 = 32;
			python.set("number1", new PyInteger(number1));
			python.set("number2", new PyInteger(number2));
			python.exec("number3 = number1+number2");
			PyObject number3 = python.get("number3");
			//System.out.println("val : "+number3.toString());
			*/
			python.execfile("practiceCode.py");
			python.exec("funcx(100)");
		}
		double endTime = System.currentTimeMillis();
		double timeSoFar1 = (endTime - startTime1)/1000;
		double timeSoFar2 = (endTime - startTime2)/1000;
		System.out.println("done");
		python.close();
		System.out.println("Time required to do complete process: "+timeSoFar1);
		System.out.println("Time required to do everything after creating the interpreter: "+timeSoFar2);
		
	}
}