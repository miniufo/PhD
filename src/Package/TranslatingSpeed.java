//
package Package;

import java.util.List;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.lagrangian.Typhoon;

//
public class TranslatingSpeed{
	//
	private static int ystr=1987;
	private static int yend=2011;
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String tranges="time=1Jan"+ystr+"-31Dec"+yend;
	
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		for(Typhoon tr:ls){
			float[] speeds=tr.getSpeeds();
			
			for(float spd:speeds) if(spd>10) System.out.println(spd);
		}
	}
}
