import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;


public class Optimize{
	public static HashMap<String, ArrayList<Circuit>> optimizer(){
		// Starts the timer
		double startTime = System.currentTimeMillis();
		double currTime = System.currentTimeMillis();
		double newStartTime = System.currentTimeMillis();
		
		//Set the costs and the operator names
		ConstantProperties.setCostPerOp();
		ConstantProperties.setOpNames();
		
		//Load values from ConstantProperties
		HashSet<String> allowedOps = ConstantProperties.allowedOps;
		double maxCost = ConstantProperties.maxCost;
		String truthValueToFind = ConstantProperties.truthValueToFind;
		String oppositeTruthValueToFind = null;
		int maxFanIn = ConstantProperties.maxFanIn;
		String dir = ConstantProperties.dir;
		int numInputs = ConstantProperties.numInputs;
		System.out.println(allowedOps);
		
		//Calculates the expected number of truth values and the size of the truth values based on the number of inputs.
		int expectedLengthOfTruthValues = 1;
		for(int i=0;i<numInputs;i++){
			expectedLengthOfTruthValues *= 2;
		}
		float expectedNumTruths = 1;
		for (int i=0;i<expectedLengthOfTruthValues;i++){
			expectedNumTruths *= 2;
		}
		System.out.println("We expect to see "+expectedNumTruths+" possible truth values");
		
		//Check if we are searching for a valid truth value.
		boolean lookingForSomething = false;
		if (truthValueToFind!=null){
			if (truthValueToFind.length()==expectedLengthOfTruthValues){
				System.out.println("Looking for "+truthValueToFind);
				lookingForSomething = true;
				oppositeTruthValueToFind = HelperFunctions.invert(truthValueToFind);
			}
			else{
				System.out.println("You entered a truth value to search for but the size of the truth value does not match the size of the truth values produced by the given number of inputs.");
				return null;
			}
			
		}
		
		//Makes sure that you included a at least the minimum number of gates to be able to find all the truth values.
		boolean isPoss = true;
		if (allowedOps.size()<4){
			isPoss = false;
			//The allowed operations must include each operator from at least one of the HashSets in combo.
			HashSet<HashSet<String>> combos = new HashSet<HashSet<String>>();
			HashSet<String> allowed = new HashSet<String>(Arrays.asList("@"));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("."));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("+","~"));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("&","~"));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList(">"));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("$"));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("+","="));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("+","^"));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("&","="));
			combos.add(allowed);
			allowed = new HashSet<String>(Arrays.asList("&","^"));
			combos.add(allowed);
			//Determines that there is a problem if it goes through each HashSet in combos and cannot find a Set for which each of
			//its operators are contained in the allowed operators.
			for (HashSet<String> acc:combos){
				boolean problem = false;
				for(String s : acc){
					if (!allowedOps.contains(s)){
						problem = true;
						break;
					}
				}
				if (!problem){
					isPoss = true;
					break;
				}
			}
		}
		//Stops the program if there it is not possible to find all truth values with the given operations.
		if (!isPoss){
			System.out.println("This combination of allowed circuits will not work.");
			return null;
		}
		
		//Adjusts the name of the output file if it is not a text file.
		Date d = new Date();
		boolean alterName = false;
		String dirNotSimple = "";
		boolean shouldDoNotSimple = false;
		if (dir.endsWith("/") || !dir.contains(".json")){
			alterName = true;
			if (!dir.endsWith("/")){
				dir += "/";
			}
			dirNotSimple = dir+d.toString().replace(' ', '_').replace(':', '_')+"FoundButNotBest.json";
			shouldDoNotSimple = true;
			String fileName = d.toString().replace(' ', '_').replace(':', '_')+"_Ops";
			
			for (String op:allowedOps){
				fileName += op;
			}
			fileName += "_MaxCost"+maxCost;
			dir = dir + fileName + ".json";
			
		}
		else if(dir.contains(".json") && dir.contains("/")){
			int indexOfSlash = dir.lastIndexOf("/");
			dirNotSimple = dir.substring(0, indexOfSlash+1)+d.toString().replace(' ', '_').replace(':', '_')+"FoundButNotBest.json";
			shouldDoNotSimple = true;
		}
		else if(dir.contains(".json") && !dir.contains("/")){
			dirNotSimple = d.toString().replace(' ', '_').replace(':', '_')+"FoundButNotBest.json";
			shouldDoNotSimple = true;
		}
		System.out.println(dir);
		
		//Make sure there are no negative costs and defaults any unspecified truth value to 1.
		//Also changes any operator with a cost of 0 that is not NOT to 1.
		ConstantProperties.setCostPerOp();
		for (String op:ConstantProperties.approvedOperators){
			if (!(ConstantProperties.costPerOp.containsKey(op))|| (ConstantProperties.costPerOp.get(op)<0) || (!op.equalsIgnoreCase("~") && ConstantProperties.costPerOp.get(op)==0)){
				ConstantProperties.costPerOp.put(op,1.0);
				System.out.println("Changing the cost of "+op+" to 1.");
			}
		}
		
		//Removes unrecognized operators
		ArrayList<String> inAllowed = new ArrayList<String> ();
		ArrayList<String> inCost = new ArrayList<String> ();
		for (String operator:ConstantProperties.costPerOp.keySet()){
			if (!ConstantProperties.approvedOperators.contains(operator)){
				inCost.add(operator);
			}
		}
		for (String operator:inCost){
			ConstantProperties.costPerOp.remove(operator);
		}
		for (String operator:allowedOps){
			if (!ConstantProperties.approvedOperators.contains(operator)){
				inAllowed.add(operator);
			}
		}
		for (String operator:inAllowed){
			allowedOps.remove(operator);
		}
		
		//Remove unused operators
		for (String operator:ConstantProperties.approvedOperators){
			if (!allowedOps.contains(operator) && ConstantProperties.costPerOp.containsKey(operator)){
				ConstantProperties.costPerOp.remove(operator);
			}
		}
		
		//Sorts the allowed operators into categories
		HashSet<String> notOps = new HashSet<String>  ();
		HashSet<String> specialZeroOps = new HashSet<String>  ();
		HashSet<String> specialOneOps = new HashSet<String>  ();
		HashSet<String> otherOps = new HashSet<String>  ();
		double costOfNot = -1;
		for (String operator:allowedOps){
			if (operator.equals("~")){
				notOps.add(operator);
				costOfNot = ConstantProperties.costPerOp.get("~");
			}
			else {
				//Every operator that is not ~ gets added here
				otherOps.add(operator);
				//These are the operators that when performed with zero invert the circuit. Only add them to the list if they will be an improvement on 
				if((operator.equals(".")||operator.equals("=")) && (!(ConstantProperties.costPerOp.containsKey("~"))||ConstantProperties.costPerOp.get(operator)<ConstantProperties.costPerOp.get("~"))){
					specialZeroOps.add(operator);
					if (ConstantProperties.costPerOp.get(operator)<costOfNot||costOfNot==-1){
						costOfNot = ConstantProperties.costPerOp.get(operator);
					}
				}
				//These are the operators that when performed with one invert the circuit. Only add them to the list if they will be an improvement on
				else if((operator.equals("@")||operator.equals("^")) && (!(ConstantProperties.costPerOp.containsKey("~"))||ConstantProperties.costPerOp.get(operator)<ConstantProperties.costPerOp.get("~"))){
					specialOneOps.add(operator);
					if (ConstantProperties.costPerOp.get(operator)<costOfNot||costOfNot==-1){
						costOfNot = ConstantProperties.costPerOp.get(operator);
					}
				}
			}
		}
		
		if (costOfNot<0){
			System.out.println("Something is wrong.");
			return null;
		}
		
		
		//These will map truth values to their minimal circuits and their cost
		HashMap<String, ArrayList<Circuit>> foundTruthValues = new HashMap<String,ArrayList<Circuit>> ();
		HashMap<String,Double> foundTruthValuesCost = new HashMap<String,Double> ();
		//These will map truth values to their minimal circuits and their cost found so far, which may not actually be the minimal.
		HashMap<String, ArrayList<Circuit>> foundTruthValuesNotSimple = new HashMap<String,ArrayList<Circuit>> ();
		HashMap<String,Double> foundTruthValuesCostNotSimple = new HashMap<String,Double> ();
		//Maps cost to a list of circuits at that cost.
		HashMap<Double, ArrayList<Circuit>> sortedByCost = new HashMap<Double,ArrayList<Circuit>> ();
		
		//Get the truth values for the one and zero circuit and for all the inputs.
		ArrayList<String> listOfTV = HelperFunctions.inputTruthValues(numInputs);
		//Create zero and one circuits
		String zeroTV = listOfTV.get(0);
		String oneTV = listOfTV.get(listOfTV.size()-1);
		Circuit zero = new Circuit("0",zeroTV);
		Circuit one = new Circuit("1",oneTV);
		ArrayList<Circuit> oneList = new ArrayList<Circuit>();
		ArrayList<Circuit> zeroList = new ArrayList<Circuit>();
		oneList.add(one);
		zeroList.add(zero);
		foundTruthValues.put(oneTV,oneList);
		foundTruthValues.put(zeroTV,zeroList);
		foundTruthValuesCost.put(oneTV,0.0);
		foundTruthValuesCost.put(zeroTV,0.0);
		foundTruthValuesNotSimple.put(oneTV,oneList);
		foundTruthValuesNotSimple.put(zeroTV,zeroList);
		foundTruthValuesCostNotSimple.put(oneTV,0.0);
		foundTruthValuesCostNotSimple.put(zeroTV,0.0);
		//Create the inputs
		ArrayList<Circuit> listOfInputs = new ArrayList<Circuit> ();
		for (int i=0;i<numInputs;i++){
			ArrayList<Circuit> tempCircList = new ArrayList<Circuit>();
			String standardInputName = "IN"+(i+1);
			String tempTV = listOfTV.get(i+1);
			Circuit tempCirc = new Circuit(standardInputName,tempTV);
			tempCircList.add(tempCirc);
			foundTruthValues.put(tempTV,tempCircList);
			foundTruthValuesCost.put(tempTV,0.0);
			foundTruthValuesNotSimple.put(tempTV,tempCircList);
			foundTruthValuesCostNotSimple.put(tempTV,0.0);
			listOfInputs.add(tempCirc);
		}
		sortedByCost.put(0.0, listOfInputs);
		
		//If there is an operation that with a cost of 0 we will run into a problem the way this program works.
		if (ConstantProperties.costPerOp.containsValue(0.0)){
			System.out.println("You cannot set the cost of something to be 0.");
			return null;
		}
		
		//initialize the cost of the first set of circuits.
		double alphaCost = 0;
		
		//This is a filtering mechanism to speed up the program. We do not want to save circuits 
		//that contain truth values that are in this list. Circuits with 0 and the minimum non-zero 
		//gate cost have their truth values added to the list.
		HashSet<String> notAllowed = new HashSet<String> ();
		int count = 0;
		 
		while(alphaCost <= maxCost && foundTruthValues.size()<expectedNumTruths){
			//Go through each cost in sortedByCost from 0 to the max. Set the group of alpha circuits.
			ArrayList<Circuit> alphaCircuits = sortedByCost.get(alphaCost);
			//Add to the list of not allowed truth values things that have a 0 smallest non-zero cost.
			if(count<=1){
				for (Circuit circ:alphaCircuits){
					notAllowed.add(circ.getTruthValue());
				}
				count++;
			}
			//Keeps track of the truth values found this round. We want to save all circuits who give these truth values (because
			//these will have the same cost as the minimum circuit found for that truth value) or whose truth values we have not
			//yet found a circuit for.
			ArrayList<String> foundThisRound = new ArrayList<String> ();
			for (Circuit circ: alphaCircuits){
				String tv = circ.getTruthValue();
				if (!foundTruthValues.containsKey(tv)){
					ArrayList<Circuit> unique = new ArrayList<Circuit> ();
					unique.add(circ);
					foundTruthValues.put(tv, unique);
					foundThisRound.add(tv);
					foundTruthValuesCost.put(tv,circ.getCost());
				}
				else if(foundTruthValues.containsKey(tv) && foundThisRound.contains(tv)){
					foundTruthValues.get(tv).add(circ);
				}
			}
			
			//Give some status updates.
			currTime = System.currentTimeMillis();
			double timeSoFar = (currTime - startTime)/1000;
			System.out.println(timeSoFar+" seconds");
			System.out.println("Number of truth values found so far: "+foundTruthValuesNotSimple.size());
			System.out.println("Number of truth values minimized so far: "+foundTruthValues.size());
			System.out.println("Potential Max Cost so far: "+Collections.max(foundTruthValuesCostNotSimple.values()));
			System.out.println("Dealing with circuits that cost "+alphaCost);
			System.out.println(alphaCircuits.size()+" circuits in this level");
			
			//Write the contents to a file.
			String tempdir = dir;
			if (alterName){
				tempdir = tempdir.replace(Double.toString(maxCost), Double.toString(alphaCost));
			}
			JsonWriter.writeToJson(foundTruthValues, foundTruthValuesCost, tempdir, timeSoFar, d, alphaCost, allowedOps, maxFanIn);
			if (shouldDoNotSimple){
				JsonWriter.writeToJson(foundTruthValuesNotSimple, foundTruthValuesCostNotSimple, dirNotSimple, timeSoFar, d, alphaCost, allowedOps, maxFanIn);
			}
			
			//Check to see if we found the truth value we were looking for if one was entered.
			if (lookingForSomething){
				boolean foundSomething = false;
				boolean foundOpposites = false;
				ArrayList<Circuit> working = new ArrayList<Circuit> ();
				ArrayList<Circuit> oppositesSameCost = new ArrayList<Circuit> ();
				ArrayList<Circuit> oppositesLess = new ArrayList<Circuit> ();
				for(String tv1 : foundTruthValues.keySet()){
					double tempCost = foundTruthValues.get(tv1).get(0).getCost();
					if (HelperFunctions.truthValuesAreSame(tv1, truthValueToFind)){
						working.addAll(foundTruthValues.get(tv1));
						foundSomething = true;
					}
					else if (HelperFunctions.truthValuesAreSame(tv1, oppositeTruthValueToFind) && tempCost==alphaCost-costOfNot){
						oppositesLess.addAll(foundTruthValues.get(tv1));
						foundOpposites = true;
					}
					else if (HelperFunctions.truthValuesAreSame(tv1, oppositeTruthValueToFind) && tempCost==alphaCost){
						oppositesSameCost.addAll(foundTruthValues.get(tv1));
						foundOpposites = true;
					}
				}
				if(foundSomething){
					System.out.println("\nFound "+truthValueToFind+":");
					System.out.println(working);
					if(foundOpposites){
						System.out.println("Also found the opposite, "+oppositeTruthValueToFind+":");
						if(oppositesLess.size()!=0){
							System.out.println(oppositesLess);
						}
						else{
							System.out.println(oppositesSameCost);
						}
					}
					System.out.println("---------------------------------------");
					return foundTruthValues;
				}
			}
			if (foundTruthValues.size() == expectedNumTruths){
				System.out.println("\nFound all truth values. Finished");
				System.out.println("---------------------------------------");
				return foundTruthValues;
			}
			if (alphaCost == maxCost){
				System.out.println("\nReached the maxCost. Finished");
				System.out.println("---------------------------------------");
				return foundTruthValues;
			}
			System.out.println("---------------------------------------");
			
			//Do all the operations that would produce the inverted circuit.
			for (Circuit circ: alphaCircuits){
				for (String tempOp : notOps){
					ArrayList<Circuit> tempArrayList = new ArrayList<Circuit> ();
					tempArrayList.add(circ);
					//Make the circuit
					Circuit tempCircuit = new Circuit(tempOp,tempArrayList);
					double tempCircuitCost = tempCircuit.getCost();
					String tempCircuitTV = tempCircuit.getTruthValue();
					if(ConstantProperties.shouldSpeedUp){
						if (foundTruthValues.containsKey(tempCircuitTV)){
							continue;
						}
					}
					//Check if this circuit is better than the one currently saved for that in foundTruthValuesNotSimple. If we have found at least one circuit for all truth values then adjust the maxCost if necessary to save space.
					if(shouldDoNotSimple && tempCircuit.canBeUsed()){
						//If this is a new truth value or the current cost for this truth value is higher than this new cost.
						if(!foundTruthValuesCostNotSimple.containsKey(tempCircuitTV)||foundTruthValuesCostNotSimple.get(tempCircuitTV)>tempCircuitCost){
							ArrayList<Circuit> tempoArrayList = new ArrayList<Circuit>();
							tempoArrayList.add(tempCircuit);
							foundTruthValuesNotSimple.put(tempCircuitTV, tempoArrayList);
							foundTruthValuesCostNotSimple.put(tempCircuitTV,tempCircuitCost);
							//Adjust maxCost if we found all 256 with less than the previous maxCost.
							if(foundTruthValuesCostNotSimple.size()==expectedNumTruths && Collections.max(foundTruthValuesCostNotSimple.values())<maxCost){
								System.out.println("A circuit has been found for all truth values, but it may not be the minimum.");
								maxCost = Collections.max(foundTruthValuesCostNotSimple.values());
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									//remove all costs that are above the new max cost
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									//trim down the circuits from the maxCost to save space.
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
						}
						//If the cost of this circuit matches the cost for the the ones we already found for the truth value
						else if(foundTruthValuesCostNotSimple.get(tempCircuitTV)==tempCircuitCost){
							foundTruthValuesNotSimple.get(tempCircuitTV).add(tempCircuit);
						}
					}
					if(lookingForSomething){
						//If we find the truth value we are looking for print it and the cost and readjust the maxCost to be the cost of that circuit.
						//Memory
						if(tempCircuit.isTruthValue(truthValueToFind) && tempCircuitCost<=maxCost){
							System.out.println("Found something with the truth value we are looking for, "+truthValueToFind);
							System.out.println(tempCircuit);
							System.out.println("Cost: "+tempCircuitCost);
							//Adjust maxCost if we do not need to go as high to find the desired truth value
							if(tempCircuitCost<maxCost){
								maxCost = tempCircuitCost;
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
						}
						//if we find the opposite of the truth value we are looking for
						//Memory
						else if(tempCircuit.isTruthValue(oppositeTruthValueToFind) && tempCircuitCost+costOfNot<=maxCost){
							System.out.println("Found something with the OPPOSITE truth value of what we are looking for, "+oppositeTruthValueToFind);
							System.out.println(tempCircuit);
							System.out.println("Cost: "+tempCircuitCost);
							System.out.println("Expected cost of actual from this: "+(tempCircuitCost+costOfNot));
							if(tempCircuitCost+costOfNot<maxCost){
								maxCost = tempCircuitCost+costOfNot;
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
						}
					}
					
					//Add circuit to sortedByCost if it is acceptable
					if(tempCircuit.canBeUsed() && !notAllowed.contains(tempCircuitTV) && tempCircuitCost<=maxCost && !(tempCircuitCost==maxCost && foundTruthValues.containsKey(tempCircuitTV))){
						//if we have circuits with that cost already, add this to the list at that cost
						 if(sortedByCost.containsKey(tempCircuitCost) && !sortedByCost.get(tempCircuitCost).contains(tempCircuit)){
							 sortedByCost.get(tempCircuitCost).add(tempCircuit);
						 }
						 //If we have not seen this cost before, make a new list with this circuit.
						 else if(!sortedByCost.containsKey(tempCircuitCost)){
							 ArrayList<Circuit> tempArray = new ArrayList<Circuit>();
							 tempArray.add(tempCircuit);
							 sortedByCost.put(tempCircuitCost,tempArray);
						 }
					}
				}
				
				for (String tempOp : specialZeroOps){
					ArrayList<Circuit> tempArrayList = new ArrayList<Circuit> ();
					tempArrayList.add(circ);
					tempArrayList.add(zero);
					Circuit tempCircuit = new Circuit(tempOp,tempArrayList);
					double tempCircuitCost = tempCircuit.getCost();
					String tempCircuitTV = tempCircuit.getTruthValue();
					if(ConstantProperties.shouldSpeedUp){
						if (foundTruthValues.containsKey(tempCircuitTV)){
							continue;
						}
					}
					if(shouldDoNotSimple&&tempCircuit.canBeUsed()){
						if(!foundTruthValuesCostNotSimple.containsKey(tempCircuitTV)||foundTruthValuesCostNotSimple.get(tempCircuitTV)>tempCircuitCost){
							ArrayList<Circuit> tempoArrayList = new ArrayList<Circuit>();
							tempoArrayList.add(tempCircuit);
							foundTruthValuesNotSimple.put(tempCircuitTV, tempoArrayList);
							foundTruthValuesCostNotSimple.put(tempCircuitTV,tempCircuitCost);
							if(foundTruthValuesCostNotSimple.size()==expectedNumTruths && Collections.max(foundTruthValuesCostNotSimple.values())<maxCost){
								System.out.println("A circuit has been found for all truth values, but it may not be the minimum.");
								maxCost = Collections.max(foundTruthValuesCostNotSimple.values());
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
							//JsonWriter.writeToJson(foundTruthValuesNotSimple, foundTruthValuesCostNotSimple, dirNotSimple, timeSoFar, d, alphaCost, allowedOps, maxFanIn);
						}
						else if(foundTruthValuesCostNotSimple.get(tempCircuitTV)==tempCircuitCost){
							foundTruthValuesNotSimple.get(tempCircuitTV).add(tempCircuit);
						}
					}
					if(lookingForSomething){
						//Memory
						if(tempCircuit.isTruthValue(truthValueToFind) && tempCircuitCost<=maxCost){
							System.out.println("Found something with the truth value we are looking for, "+truthValueToFind);
							System.out.println(tempCircuit);
							System.out.println("Cost: "+tempCircuitCost);
							if(tempCircuitCost<maxCost){
								maxCost = tempCircuitCost;
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
						}
						//Memory
						else if(tempCircuit.isTruthValue(oppositeTruthValueToFind) && tempCircuitCost+costOfNot<=maxCost){
							System.out.println("Found something with the OPPOSITE truth value of what we are looking for, "+oppositeTruthValueToFind);
							System.out.println(tempCircuit);
							System.out.println("Cost: "+tempCircuitCost);
							System.out.println("Expected cost of actual from this: "+(tempCircuitCost+costOfNot));
							if(tempCircuitCost+costOfNot<maxCost){
								maxCost = tempCircuitCost+costOfNot;
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
						}
					}
					
					//Add circuit to sortedByCost if it is acceptable
					if(tempCircuit.canBeUsed() && !notAllowed.contains(tempCircuitTV) && tempCircuitCost<=maxCost && !(tempCircuitCost==maxCost && foundTruthValues.containsKey(tempCircuitTV))){
						//if we have circuits with that cost already, add this to the list at that cost
						 if(sortedByCost.containsKey(tempCircuitCost) && !sortedByCost.get(tempCircuitCost).contains(tempCircuit)){
							 sortedByCost.get(tempCircuitCost).add(tempCircuit);
						 }
						 //If we have not seen this cost before, make a new list with this circuit.
						 else if(!sortedByCost.containsKey(tempCircuitCost)){
							 ArrayList<Circuit> tempArray = new ArrayList<Circuit>();
							 tempArray.add(tempCircuit);
							 sortedByCost.put(tempCircuitCost,tempArray);
						 }
					}
				}
				
				for (String tempOp : specialOneOps){
					ArrayList<Circuit> tempArrayList = new ArrayList<Circuit> ();
					tempArrayList.add(circ);
					tempArrayList.add(one);
					Circuit tempCircuit = new Circuit(tempOp,tempArrayList);
					double tempCircuitCost = tempCircuit.getCost();
					String tempCircuitTV = tempCircuit.getTruthValue();
					if(ConstantProperties.shouldSpeedUp){
						if (foundTruthValues.containsKey(tempCircuitTV)){
							continue;
						}
					}
					if(shouldDoNotSimple&&tempCircuit.canBeUsed()){
						if(!foundTruthValuesCostNotSimple.containsKey(tempCircuitTV)||foundTruthValuesCostNotSimple.get(tempCircuitTV)>tempCircuitCost){
							ArrayList<Circuit> tempoArrayList = new ArrayList<Circuit>();
							tempoArrayList.add(tempCircuit);
							foundTruthValuesNotSimple.put(tempCircuitTV, tempoArrayList);
							foundTruthValuesCostNotSimple.put(tempCircuitTV,tempCircuitCost);
							if(foundTruthValuesCostNotSimple.size()==expectedNumTruths && Collections.max(foundTruthValuesCostNotSimple.values())<maxCost){
								System.out.println("A circuit has been found for all truth values, but it may not be the minimum.");
								maxCost = Collections.max(foundTruthValuesCostNotSimple.values());
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
							//JsonWriter.writeToJson(foundTruthValuesNotSimple, foundTruthValuesCostNotSimple, dirNotSimple, timeSoFar, d, alphaCost, allowedOps, maxFanIn);
						}
						else if(foundTruthValuesCostNotSimple.get(tempCircuitTV)==tempCircuitCost){
							foundTruthValuesNotSimple.get(tempCircuitTV).add(tempCircuit);
						}
					}
					if(lookingForSomething){
						//Memory
						if(tempCircuit.isTruthValue(truthValueToFind) && tempCircuitCost<=maxCost){
							System.out.println("Found something with the truth value we are looking for, "+truthValueToFind);
							System.out.println(tempCircuit);
							System.out.println("Cost: "+tempCircuitCost);
							if(tempCircuitCost<maxCost){
								maxCost = tempCircuitCost;
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
						}
						//Memory
						else if(tempCircuit.isTruthValue(oppositeTruthValueToFind) && tempCircuitCost+costOfNot<=maxCost){
							System.out.println("Found something with the OPPOSITE truth value of what we are looking for, "+oppositeTruthValueToFind);
							System.out.println(tempCircuit);
							System.out.println("Cost: "+tempCircuitCost);
							System.out.println("Expected cost of actual from this: "+(tempCircuitCost+costOfNot));
							if(tempCircuitCost+costOfNot<maxCost){
								maxCost = tempCircuitCost+costOfNot;
								System.out.println("Changing maxCost to "+maxCost);
								Set<Double> allCostsFound = sortedByCost.keySet();
								ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
								allCostsFoundDup.addAll(allCostsFound);
								for(Double tcost : allCostsFoundDup){
									if(tcost>maxCost && sortedByCost.containsKey(tcost)){
										sortedByCost.remove(tcost);
									}
									else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
										ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
										tempCopyList.addAll(sortedByCost.get(tcost));
										for (Circuit tCircuit : tempCopyList){
											if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
												sortedByCost.get(tcost).remove(tCircuit);
											}
										}
									}
								}
							}
						}
					}
					
					//Add circuit to sortedByCost if it is acceptable
					if(tempCircuit.canBeUsed() && !notAllowed.contains(tempCircuitTV) && tempCircuitCost<=maxCost && !(tempCircuitCost==maxCost && foundTruthValues.containsKey(tempCircuitTV))){
						//if we have circuits with that cost already, add this to the list at that cost
						 if(sortedByCost.containsKey(tempCircuitCost) && !sortedByCost.get(tempCircuitCost).contains(tempCircuit)){
							 sortedByCost.get(tempCircuitCost).add(tempCircuit);
						 }
						 //If we have not seen this cost before, make a new list with this circuit.
						 else if(!sortedByCost.containsKey(tempCircuitCost)){
							 ArrayList<Circuit> tempArray = new ArrayList<Circuit>();
							 tempArray.add(tempCircuit);
							 sortedByCost.put(tempCircuitCost,tempArray);
						 }
					}
				}
			}
			
			
			
			//Get every unique combination of circuits with a number of circuits equal to or less than the maxFanIn.
			ArrayList<ArrayList<Circuit>> allCombinationsOfCircuits = new ArrayList<ArrayList<Circuit>>();
			ArrayList<Circuit> allCircuitsCombined = new ArrayList<Circuit> ();
			for (Double tcost : sortedByCost.keySet()){
				if (tcost<=alphaCost){
					allCircuitsCombined.addAll(sortedByCost.get(tcost));
				}
			}
			for (int i=1;i<maxFanIn;i++){
				allCombinationsOfCircuits.addAll((Collection<? extends ArrayList<Circuit>>) HelperFunctions.combinations(allCircuitsCombined, i, false));
			}
			//Make the circuits from the combinations made
			for (Circuit circ : alphaCircuits){
				for(ArrayList<Circuit> tempArrayOfCircuits : allCombinationsOfCircuits){
					//Print progress every 5 minutes
					if(((System.currentTimeMillis() - newStartTime)/60000)>=5){
						newStartTime = System.currentTimeMillis();
						double totalProgress = ((100.0*alphaCircuits.indexOf(circ))/alphaCircuits.size()) + (100.0*(allCombinationsOfCircuits.indexOf(tempArrayOfCircuits)+1)/alphaCircuits.size())/allCombinationsOfCircuits.size();
						currTime = System.currentTimeMillis();
						timeSoFar = (currTime - startTime)/60000;
						JsonWriter.writeToJson(foundTruthValuesNotSimple, foundTruthValuesCostNotSimple, dirNotSimple, timeSoFar, d, alphaCost, allowedOps, maxFanIn);
						System.out.println(totalProgress + "% complete with building from "+ alphaCost+" cost gates after "+timeSoFar+" minutes");
						System.out.println("Number of truth values found so far: "+foundTruthValuesNotSimple.size());
						System.out.println("Number of truth values minimized so far: "+foundTruthValues.size());
						System.out.println("Potential Max Cost so far: "+Collections.max(foundTruthValuesCostNotSimple.values()));
					}
					//We do not want to perform any operations with a group of circuits if it contains the alpha circuit because we will just get 1,0, or the inverse of the alpha circuit and we already took care of all the inverses.
					if(tempArrayOfCircuits.contains(circ)){
						continue;
					}
					//Make a copy so we aren't changing the original
					ArrayList<Circuit> tempArrayOfCircuitsCopy = new ArrayList<Circuit> ();
					tempArrayOfCircuitsCopy.add(circ);
					tempArrayOfCircuitsCopy.addAll(tempArrayOfCircuits);
					for(String tempOp : otherOps){
						Circuit tempCircuit = new Circuit(tempOp,tempArrayOfCircuitsCopy);
						double tempCircuitCost = tempCircuit.getCost();
						String tempCircuitTV = tempCircuit.getTruthValue();
						if(ConstantProperties.shouldSpeedUp){
							if (foundTruthValues.containsKey(tempCircuitTV)){
								continue;
							}
						}
						if(shouldDoNotSimple && tempCircuit.canBeUsed()){
							if(!foundTruthValuesCostNotSimple.containsKey(tempCircuitTV)||foundTruthValuesCostNotSimple.get(tempCircuitTV)>tempCircuitCost){
								ArrayList<Circuit> tempoArrayList = new ArrayList<Circuit>();
								tempoArrayList.add(tempCircuit);
								foundTruthValuesNotSimple.put(tempCircuitTV, tempoArrayList);
								foundTruthValuesCostNotSimple.put(tempCircuitTV,tempCircuitCost);
								if(foundTruthValuesCostNotSimple.size()==expectedNumTruths && Collections.max(foundTruthValuesCostNotSimple.values())<maxCost){
									System.out.println("A circuit has been found for all truth values, but it may not be the minimum.");
									maxCost = Collections.max(foundTruthValuesCostNotSimple.values());
									System.out.println("Changing maxCost to "+maxCost);
									Set<Double> allCostsFound = sortedByCost.keySet();
									ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
									allCostsFoundDup.addAll(allCostsFound);
									for(Double tcost : allCostsFoundDup){
										if(tcost>maxCost && sortedByCost.containsKey(tcost)){
											sortedByCost.remove(tcost);
										}
										else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
											ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
											tempCopyList.addAll(sortedByCost.get(tcost));
											for (Circuit tCircuit : tempCopyList){
												if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
													sortedByCost.get(tcost).remove(tCircuit);
												}
											}
										}
									}
								}
								//JsonWriter.writeToJson(foundTruthValuesNotSimple, foundTruthValuesCostNotSimple, dirNotSimple, timeSoFar, d, alphaCost, allowedOps, maxFanIn);
							}
							else if(foundTruthValuesCostNotSimple.get(tempCircuitTV)==tempCircuitCost){
								foundTruthValuesNotSimple.get(tempCircuitTV).add(tempCircuit);
							}
						}
						if(lookingForSomething){
							//Memory
							if(tempCircuit.isTruthValue(truthValueToFind) && tempCircuitCost<=maxCost){
								System.out.println("Found something with the truth value we are looking for, "+truthValueToFind);
								System.out.println(tempCircuit);
								System.out.println("Cost: "+tempCircuitCost);
								if(tempCircuitCost<maxCost){
									maxCost = tempCircuitCost;
									System.out.println("Changing maxCost to "+maxCost);
									Set<Double> allCostsFound = sortedByCost.keySet();
									ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
									allCostsFoundDup.addAll(allCostsFound);
									for(Double tcost : allCostsFoundDup){
										if(tcost>maxCost && sortedByCost.containsKey(tcost)){
											sortedByCost.remove(tcost);
										}
										else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
											ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
											tempCopyList.addAll(sortedByCost.get(tcost));
											for (Circuit tCircuit : tempCopyList){
												if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
													sortedByCost.get(tcost).remove(tCircuit);
												}
											}
										}
									}
								}
							}
							//Memory
							else if(tempCircuit.isTruthValue(oppositeTruthValueToFind) && tempCircuitCost+costOfNot<=maxCost){
								System.out.println("Found something with the OPPOSITE truth value of what we are looking for, "+oppositeTruthValueToFind);
								System.out.println(tempCircuit);
								System.out.println("Cost: "+tempCircuitCost);
								System.out.println("Expected cost of actual from this: "+(tempCircuitCost+costOfNot));
								if(tempCircuitCost+costOfNot<maxCost){
									maxCost = tempCircuitCost+costOfNot;
									System.out.println("Changing maxCost to "+maxCost);
									Set<Double> allCostsFound = sortedByCost.keySet();
									ArrayList<Double> allCostsFoundDup = new ArrayList<Double> ();
									allCostsFoundDup.addAll(allCostsFound);
									for(Double tcost : allCostsFoundDup){
										if(tcost>maxCost && sortedByCost.containsKey(tcost)){
											sortedByCost.remove(tcost);
										}
										else if(tcost==maxCost && sortedByCost.containsKey(tcost)){
											ArrayList<Circuit> tempCopyList = new ArrayList<Circuit> ();
											tempCopyList.addAll(sortedByCost.get(tcost));
											for (Circuit tCircuit : tempCopyList){
												if(foundTruthValues.containsKey(tCircuit.getTruthValue())){
													sortedByCost.get(tcost).remove(tCircuit);
												}
											}
										}
									}
								}
							}
						}
						
						//Add circuit to sortedByCost if it is acceptable
						if(tempCircuit.canBeUsed() && !notAllowed.contains(tempCircuitTV) && tempCircuitCost<=maxCost && !(tempCircuitCost==maxCost && foundTruthValues.containsKey(tempCircuitTV))){
							//if we have circuits with that cost already, add this to the list at that cost
							 if(sortedByCost.containsKey(tempCircuitCost) && !sortedByCost.get(tempCircuitCost).contains(tempCircuit)){
								 sortedByCost.get(tempCircuitCost).add(tempCircuit);
							 }
							 //If we have not seen this cost before, make a new list with this circuit.
							 else if(!sortedByCost.containsKey(tempCircuitCost)){
								 ArrayList<Circuit> tempArray = new ArrayList<Circuit>();
								 tempArray.add(tempCircuit);
								 sortedByCost.put(tempCircuitCost,tempArray);
							 }
						}
					}
				}
			}
			
			//Gets the next cost for the alpha circuits.
			ArrayList<Double> allCosts =new ArrayList<Double> (sortedByCost.keySet());
			Collections.sort(allCosts);
			//If we reach the last cost available, break.
			if(allCosts.indexOf(alphaCost)==allCosts.size()-1){
				break;
			}
			alphaCost = allCosts.get(allCosts.indexOf(alphaCost)+1);
		}
		System.out.println("Finished");
		return foundTruthValues;
	}
	
}