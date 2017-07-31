package Package;
//
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;


//
public class CylindricStations{
	//
	private static final float threshold=30f;
	
	private static final boolean noDepress=true;	// wind > 17.2 m/s
	private static final boolean noLanding=true;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String respath="d:/Data/PhD/Climatology/"+dsets+"/";
	private static final String tranges="time=1Jan1989-31Dec2009";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		DiagnosisFactory df2=DiagnosisFactory.parseFile("d:/Data/ERAInterim/lsm.ctl");
		
		CtlDescriptor ctl2=(CtlDescriptor)df2.getDataDescriptor();
		Variable lsm=df2.getVariables(new Range("z(1,1)",ctl2),false,"lsm")[0];
		
		int cc=0;
		for(Typhoon tr:ls){
			if(++cc%20==0) System.out.print(".");
			
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString(36,25,2,0.3f,-650));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v");
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getZonalVelocity(),tr.getMeridionalVelocity(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable efcsm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),9,18);	// 300-600 km
			
			Variable flux=utvr[0].copy();	flux.multiplyEq(utvr[1]);
			flux.setName("flx");	flux.setCommentAndUnit("eddy momentum flux, u'v' (m^2 s^-2)");
			
			ct=new CoordinateTransformation(new SphericalSpatialModel(ctl2),csm);
			Variable lsmm=dm.cRadialAverage(ct.transToCylindricalInvariantly(lsm),1,6).anomalizeX();
			
			boolean[] wind=noDepress?	// wind >= 17.2 m/s
				IntensityModel.greaterEqualThan(tr.getWinds(),17.2f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] lsmb=noLanding?	// no land within 200 km
				IntensityModel.lessThan(lsmm.getData()[0][0][0],1e-9f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] llon=withinWNP?	// lons >= 100E
				IntensityModel.greaterEqualThan(tr.getLongitudes(),100):IntensityModel.newBooleans(tr.getTCount());
			boolean[] rlon=withinWNP?	// lons <= 190E
				IntensityModel.lessEqualThan(tr.getLongitudes(),190):IntensityModel.newBooleans(tr.getTCount());
			
			boolean[] hghEFCS =IntensityModel.ABSgreaterEqualThan(efcsm.getData()[1][0][0],threshold);
			boolean[] combEFCS=IntensityModel.combination(wind,hghEFCS,lsmb,llon,rlon);
			
			if(IntensityModel.getValidCount(combEFCS)>0){
				int year=(int)(tr.getTime(0)/10000000000L);
				int month=(int)(tr.getTime(0)%10000000000L/100000000L);
				
				String filename=year+String.format("%02d",month)+tr.getName()+".dat";
				String ctlname=year+String.format("%02d",month)+tr.getName()+".ctl";
				
				DataWrite dw=DataIOFactory.getDataWrite(dd,respath+"station/"+filename);
				dw.writeData(dd,utvr[0],utvr[1],flux);	dw.closeFile();
				
				try{
					BufferedReader br=new BufferedReader(
						new InputStreamReader(
							Runtime.getRuntime().exec("stnmap -i "+respath+"station/"+ctlname).getInputStream()
						)
					);
					
					boolean print=false;
					
					if(print){
						String str=null;
						
						System.out.println();
						while((str=br.readLine())!=null) System.out.println(str);
						System.out.println();
						
						br.close();
					}
					
				}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			}
		}
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("\nthreshold: "+threshold+" m s^-1 / day");
	}
}
