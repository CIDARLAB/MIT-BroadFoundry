import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class ConstantProperties{
	//Do not change these.
	public static ArrayList<String> approvedOperators = new ArrayList<String> (Arrays.asList("~", "&", "@","+","^",".","="));
	public static HashMap<String,Double> costPerOp = new HashMap<String,Double>();
	public static HashMap<String,String> opNames = new HashMap<String,String> ();
	
	//The number of inputs can only be two or three.
	public static int numInputs = 3;
	//Choose which operators you want to include.
	public static HashSet<String> allowedOps= new HashSet<String>(Arrays.asList(".","~"));
	//Choose a max cost you want to go up to.
	public static double maxCost = 8;
	//Choose the maxFanIn for each gate;
	public static int maxFanIn = 2;
	//If you are searching for a specific truth value, enter it here. Otherwise leave it as null.
	public static String truthValueToFind = null;
	//Enter the location of where you want to write the Json file.
	public static String dir = "C:/Users/Arinze/OneDrive/UROP_Summer_2015/NewOptimizer";
	//Set this to true if you want slightly less accurate results that are obtained slightly faster
	public static boolean shouldSpeedUp = false;

	//Set the cost for each operator as you see fit. Only NOT can have a cost of 0. All must be positive
	public static void setCostPerOp(){
		costPerOp.put("~",1.0);
		costPerOp.put("&",1.0);
		costPerOp.put("@",1.0);
		costPerOp.put("+",1.0);
		costPerOp.put("^",1.0);
		costPerOp.put(".",1.0);
		costPerOp.put("=",1.0);
	}
	
	public static void setOpNames(){
		opNames.put("~","NOT");
		opNames.put("&","AND");
		opNames.put("@","NAND");
		opNames.put("+","OR");
		opNames.put("^","XOR");
		opNames.put(".","NOR");
		opNames.put("=","XNOR");
	}
}