package Package;

import miniufo.database.AccessBestTrack.DataSets;
import static java.lang.Math.exp;

//
public final class IntensityModel{
	//
	public static final float undef=-9999;
	
	
	// compute maximum potential intensity based on SST
	public static float[] cMPI(float[] sst){
		float[] re=sst.clone();
		
		for(int i=0,I=sst.length;i<I;i++)
		re[i]=15.69f+(float)(98.03*exp(0.1806*(sst[i]-30.0)));
		
		return re;
	}
	
	// compute potential intensity
	public static float[] cPotential(float[] sst,float[] wnd){
		float[] re=sst.clone();
		
		for(int i=0,I=sst.length;i<I;i++){
			re[i]=15.69f+(float)(98.03*exp(0.1806*(sst[i]-30.0))); // Zeng et al. 2007 MWR
			//re[i]=(float)(74.0*exp(0.2*(sst[i]-25.0))); // DeMaria et al. 1993 JAS
			re[i]-=wnd[i]==0?15:wnd[i];
		}
		
		return re;
	}
	
	// compute the atan according to the wind (degree)
	public static float[] cAtanByWind(float[] u,float[] v){
		int t=u.length;
		
		float[] re=new float[t];
		
		for(int l=0;l<t;l++)
		re[l]=(float)Math.toDegrees(Math.atan2(v[l],u[l]));
		
		return re;
	}
	
	// count of the valid records
	public static int getValidCount(boolean[] bs){
		int count=0;
		
		for(boolean b:bs)
		if(b) count++;
		
		return count;
	}
	
	//
	public static int getCycles(boolean[] data){
		int count=0; boolean tmp=false;
		
		if(data[0]) return -1;	// TS first
		
		for(boolean b:data)
		if(b!=tmp){
			count++;
			tmp=b;
		}
		
		if(count%2==1){
			if(count>1) return -2;	// odd count
			else count++;
		}
		
		return count/2;
	}
	
	
	// compile boolean arrays
	public static boolean[] combination(boolean[]... bs){
		if(bs==null||bs.length<1) throw new IllegalArgumentException("invalid args");
		
		int N=bs.length,T=bs[0].length;
		
		for(int i=1;i<N;i++) if(T!=bs[i].length)
		throw new IllegalArgumentException("lengths not equal");
		
		boolean[] re=bs[0].clone();
		
		for(int l=0;l<T;l++)
		for(int i=0;i<N;i++)
		re[l]&=bs[i][l];
		
		return re;
	}
	
	
	// new boolean[] with default true
	public static boolean[] newBooleans(int n){
		boolean[] re=new boolean[n];
		
		for(int i=0;i<n;i++) re[i]=true;
		
		return re;
	}
	
	
	// whether the angles is in the range (degree)
	public static boolean[] inRange(float[] angles,float center,float radius){
		int n=angles.length;
		
		boolean[] re=new boolean[n];
		
		if(center-radius>-180&&center+radius<180){
			for(int l=0;l<n;l++)
			if((center-radius)<=angles[l]&&angles[l]<(center+radius))
			re[l]=true;
			
		}else{
			for(int l=0;l<n;l++)
			if(Math.abs(angles[l])>=135+22.5f)
			re[l]=true;
		}
		
		return re;
	}
	
	
	//whether data is within threshold
	public static boolean[] greaterEqualThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(data[i]>=threshold) del[i]=true;
		
		return del;
	}
	
	public static boolean[] lessEqualThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(data[i]<=threshold) del[i]=true;
		
		return del;
	}
	
	public static boolean[] greaterThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(data[i]>threshold) del[i]=true;
		
		return del;
	}
	
	public static boolean[] lessThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(data[i]<threshold) del[i]=true;
		
		return del;
	}
	
	public static boolean[] ABSgreaterEqualThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(Math.abs(data[i])>=threshold) del[i]=true;
		
		return del;
	}
	
	public static boolean[] ABSgreaterThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(Math.abs(data[i])>threshold) del[i]=true;
		
		return del;
	}
	
	public static boolean[] ABSlessEqualThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(Math.abs(data[i])<=threshold) del[i]=true;
		
		return del;
	}
	
	public static boolean[] ABSlessThan(float[] data,float threshold){
		int N=data.length;
		
		boolean[] del=new boolean[N];
		
		for(int i=0;i<N;i++)
		if(Math.abs(data[i])<threshold) del[i]=true;
		
		return del;
	}
	
	
	/**
	 * Compute changes of the data (delta-data) using centered time differencing.
	 * The first and last are computed using one-sided time differencing.
	 */
	public static float[] getChangesByCentralDiff(float[] data){
		int N=data.length;
		
		if(N==1) return new float[1];
		
		float[] ch=new float[N];
		
		ch[0]=data[1]-data[0];
		
		for(int l=1,L=N-1;l<L;l++) ch[l]=data[l+1]-data[l-1];
		
		ch[N-1]=data[N-1]-data[N-2];
		
		return ch;
	}
	
	/**
	 * Compute changes of the data (delta-data) using forward time differencing.
	 */
	public static float[] getChangesByForwardDiff(float[] data,int interv){
		int N=data.length;
		
		float[] ch=new float[N];
		
		for(int l=0;l<N;l++) ch[l]=undef;
		
		for(int l=0,L=N-interv;l<L;l++) ch[l]=data[l+interv]-data[l];
		
		return ch;
	}
	
	
	/**
     * get TC data
     */
	public static String getPath(DataSets ds){
		final String CMAPath ="D:/Data/Typhoons/CMA/CMA.txt";
		final String JMAPath ="D:/Data/Typhoons/JMA/JMA.txt";
		final String JTWCPath="D:/Data/Typhoons/JTWC/JTWC.txt";
		
		switch(ds){
			case CMA:  return CMAPath;
			case JMA:  return JMAPath;
			case JTWC: return JTWCPath;
			default:   throw new IllegalArgumentException("unsupported DataSet");
		}
	}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(Math.toDegrees(Math.atan2(0,1)));	// ww
		System.out.println(Math.toDegrees(Math.atan2(1,1)));	// nw
		System.out.println(Math.toDegrees(Math.atan2(1,1)));	// nn
		System.out.println(Math.toDegrees(Math.atan2(1,1)));	// ne
		System.out.println(Math.toDegrees(Math.atan2(1,1)));	// ee
		System.out.println(Math.toDegrees(Math.atan2(1,1)));	// se
		System.out.println(Math.toDegrees(Math.atan2(1,1)));	// ss
		System.out.println(Math.toDegrees(Math.atan2(1,1)));	// sw
	}*/
}
