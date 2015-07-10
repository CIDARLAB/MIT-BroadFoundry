public class FindMinimalCircuits{
	public static void main(String [] args){

		Optimize.optimizer1();
		ConstantProperties.allowedOps.remove(".");
		ConstantProperties.allowedOps.add("@");
		Optimize.optimizer1();
		ConstantProperties.allowedOps.remove("@");
		ConstantProperties.allowedOps.add("+");
		ConstantProperties.allowedOps.add("&");
		ConstantProperties.allowedOps.add("~");
		//ConstantProperties.maxCost = 8;
		Optimize.optimizer1();
	}
	

}