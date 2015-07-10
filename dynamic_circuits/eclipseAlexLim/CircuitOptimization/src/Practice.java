import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;

import javax.rmi.CORBA.Util;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class Practice{
	public static void main(String[] args){
		/*HashMap<String,ArrayList<Circuit>> answer = Optimize.optimizer1();
		ArrayList<JSONObject> x = new ArrayList<JSONObject>();
		JSONObject y = new JSONObject();
		ArrayList<String> truthValues = new ArrayList<String>(answer.keySet());
		Collections.sort(truthValues);
		for (String tv:truthValues){
			y = new JSONObject();
			JSONArray temp = new JSONArray();
			ArrayList<Circuit> temp1 = answer.get(tv);
			ArrayList<String> temp2 = new ArrayList<String>();
			for(Circuit i:temp1){
				temp2.add(i.toString());
			}
			temp.addAll(temp2);
			y.put("truthValue", tv);
			y.put("circuits", temp);
			x.add(y);
			
		}

		String filepath = "C:/Users/Arinze/SkyDrive/Eclipse Workspace/stuff.json";
		prettyPrintJSONArray(x,filepath, false);
		
		Date d = new Date();
		System.out.println(d);
		String x = d.toString().replace(' ', '_');
		System.out.println(x);*/
	}
	public static void prettyPrintJSONArray(ArrayList<JSONObject> objects, String filepath, boolean append) {


        String json_array = "";


        for(int i=0; i<objects.size(); ++i) {


            JSONObject obj = objects.get(i);


            Gson gson = new GsonBuilder().disableHtmlEscaping().setPrettyPrinting().create();
            String json = gson.toJson(obj);


            json_array += json;



            if(i < objects.size() - 1) {
                json_array += ",\n";
            }



        }

        JsonWriter.fileWriter(filepath, "[\n" + json_array + "\n]\n", append);
    }
	
}