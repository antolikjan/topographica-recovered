import __main__
import numpy
import pylab
import scipy.stats
import os.path
import os
import copy
import pickle
import param
from math import pi, sqrt, exp, pow
from topo.plotting.plotfilesaver import * 
from topo.command.pylabplot import cyclic_tuning_curve, matrixplot
from topo.command.analysis import save_plotgroup
from matplotlib.ticker import MaxNLocator
from param import normalize_path
from topo.analysis.featureresponses import MeasureResponseCommand, FeatureMaps, FeatureCurveCommand, UnitCurveCommand, SinusoidalMeasureResponseCommand,PatternPresenter, FeatureResponses
import matplotlib.gridspec as gridspec
from matplotlib import rc
import topo
from numpy.random import RandomState

def circular_dist(a, b, period):
    """
    Returns the distance between a and b (scalars) in a domain with `period` period.
    """
    return  numpy.minimum(numpy.abs(a - b), period - numpy.abs(a - b))


#rc('text', usetex=True)
rc('mathtext',default='regular')
pylab.rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=8)    
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
    'weight' : 'normal', 'size' : 0}
rc('xtick', labelsize=20)
rc('ytick', labelsize=20)
rc('legend',fontsize=20)
rc('legend',labelspacing=0.25)

def disable_top_right_axis(ax):
    for loc, spine in ax.spines.iteritems():
            if loc in ['right','top']:
               spine.set_color('none') # don't draw spine
    for tick in ax.yaxis.get_major_ticks():
        tick.tick2On = False
    for tick in ax.yaxis.get_minor_ticks():
        tick.tick2On = False
    for tick in ax.xaxis.get_major_ticks():
        tick.tick2On = False
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick2On = False

def disable_bottom_axis(ax):
    for loc, spine in ax.spines.iteritems():
            if loc in ['bottom']:
               spine.set_color('none') # don't draw spine
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1On = False

def disable_left_axis(ax):
    for loc, spine in ax.spines.iteritems():
            if loc in ['left']:
               spine.set_color('none') # don't draw spine
    for tick in ax.yaxis.get_major_ticks():
        tick.tick1On = False

def remove_x_tick_labels():
    pylab.xticks([],[])  
 
def remove_y_tick_labels(): 
    pylab.yticks([],[])  

class surround_analysis():

    peak_near_facilitation_hist = []
    peak_supression_hist  = []   
    peak_far_facilitation_hist  = []
    sheet_name = ""
    data_dict = {}
    
    low_contrast=__main__.__dict__.get('LC',110)
    high_contrast=__main__.__dict__.get('HC',230)
    
    def __init__(self,sheet_name="V1Complex",prefix=None):
        if prefix:
            f = open(prefix+'data_dict.pickle','rb')
            (self.OR,self.OS,self.MR,self.data_dict) = pickle.load(f)
            f.close()
            
            if True:
                self.lhi = compute_local_homogeneity_index(self.OR*pi,__main__.__dict__.get('LHI',2.0))    
                f = open(normalize_path('lhi'+str(__main__.__dict__.get('LHI',2.0))+'.pickle'),'wb')            
                pickle.dump(self.lhi,f)
                f.close()
            else:        
                f = open(normalize_path('lhi'+str(__main__.__dict__.get('LHI',2.0))+'.pickle'),'rb')            
                self.lhi = pickle.load(f)
        else:
            import topo
            self.sheet_name=sheet_name
            self.sheet=topo.sim[sheet_name]
            # Center mask to matrixidx center
            self.center_r,self.center_c = self.sheet.sheet2matrixidx(0,0)
            self.center_x,self.center_y = self.sheet.matrixidx2sheet(self.center_r,self.center_c)
        
        
            from topo.analysis.featureresponses import PatternPresenter            
            PatternPresenter.duration=__main__.__dict__.get('duration',4.0) #!
            import topo.command.pylabplot
            reload(topo.command.pylabplot)

        
            FeatureCurveCommand.display=True
            FeatureCurveCommand.curve_parameters=[{"contrast":self.high_contrast}]
            FeatureCurveCommand.sheet=topo.sim[sheet_name]
            SinusoidalMeasureResponseCommand.num_phase=8
            SinusoidalMeasureResponseCommand.frequencies=[__main__.__dict__.get('FREQ',2.4)]
            SinusoidalMeasureResponseCommand.scale=self.high_contrast/100.0
            MeasureResponseCommand.scale=self.high_contrast/100.0
            FeatureCurveCommand.num_orientation=__main__.__dict__.get('num_orientations',12)
            FeatureResponses.repetitions = __main__.__dict__.get('repetitions',1)
            
	    import topo.analysis.vision
            topo.analysis.vision.measure_and_analyze_complexity.instance(sheet=topo.sim["V1Complex"])()
	    save_plotgroup("Orientation Preference and Complexity")
            #topo.command.pylabplot.measure_or_tuning_fullfield.instance(sheet=topo.sim["V1Complex"])()
            
            FeatureCurveCommand.curve_parameters=[{"contrast":self.low_contrast},{"contrast":self.high_contrast}]
            
            self.OR = topo.sim["V1Complex"].sheet_views['OrientationPreference'].view()[0]
            self.OS = topo.sim["V1Complex"].sheet_views['OrientationSelectivity'].view()[0]
            self.MR = topo.sim["V1Complex"].sheet_views['ComplexSelectivity'].view()[0]

            
       
    def run_analysis_with_step_grid(self,grid_step_radius,step_size,max_curves=None):
        steps = []
        for i in xrange(0,grid_step_radius*2+1):
            for j in xrange(0,grid_step_radius*2+1):
                steps.append([self.center_r+(i-grid_step_radius)*step_size,self.center_c+(j-grid_step_radius)*step_size])
        if max_curves != None:
            self.analyse(steps[0:max_curves],ns=__main__.__dict__.get('number_sizes',10))
        else:
            self.analyse(steps,ns=__main__.__dict__.get('number_sizes',10))
            
    def run_lhi_informed_analysis(self,max_curves=26,center_size=20,index=None):
	print self.OR
        if True:
            self.lhi = compute_local_homogeneity_index(self.OR*pi,__main__.__dict__.get('LHI',2.0))
            f = open(normalize_path('lhi'+str(__main__.__dict__.get('LHI',2.0))+'.pickle'),'wb')
            pickle.dump(self.lhi,f)
            f.close()
        else:
            f = open(normalize_path('lhi'+str(__main__.__dict__.get('LHI',2.0))+'.pickle'),'rb')            
            self.lhi = pickle.load(f)
        
        lhi_center = self.lhi[self.center_r-center_size:self.center_r+center_size,self.center_c-center_size:self.center_c+center_size]
	steps = []
        
        r = RandomState(1023)
       
	if not __main__.__dict__.get('uniform',False):

 
	        pinwheels = r.permutation(numpy.nonzero(numpy.ravel(lhi_center) < __main__.__dict__.get('cutoff',0.3))[0])
        	domains = r.permutation(numpy.nonzero(numpy.ravel(lhi_center) > (1-__main__.__dict__.get('cutoff',0.3)))[0])
	        
		assert len(pinwheels) > max_curves/2

        	#s = numpy.argsort(numpy.ravel(lhi_center))
        
	        if index == None:
	           for i in xrange(0,max_curves/2):
        	        (x,y) = numpy.unravel_index(pinwheels[i],lhi_center.shape)
                	steps.append((x+self.center_r-center_size,y+self.center_c-center_size))

	                (x,y) = numpy.unravel_index(domains[i],lhi_center.shape)
        	        steps.append((x+self.center_r-center_size,y+self.center_c-center_size))
	        else:
                	if (index % 2) == 0:
        	           (x,y) = numpy.unravel_index(pinwheels[int(index/2)],lhi_center.shape)
	                   steps= [(x+self.center_r-center_size,y+self.center_c-center_size)]
        	        else:
                	   (x,y) = numpy.unravel_index(domains[int(index/2)],lhi_center.shape)
	                   steps= [(x+self.center_r-center_size,y+self.center_c-center_size)] 
        else:
		bins = []
		for i in xrange(0,10):
		    a = numpy.ravel(lhi_center) >= i*0.1	
		    b = numpy.ravel(lhi_center) <  (i+1)*0.1
		    bins.append(r.permutation(numpy.nonzero(numpy.multiply(a,b))[0]))
		
		(x,y) = numpy.unravel_index(bins[index % 10][int(index/10)],lhi_center.shape)
                steps= [(x+self.center_r-center_size,y+self.center_c-center_size)]
   
	    	    					

		#places = r.permutation(numpy.arange(0,len(numpy.ravel(lhi_center)),1))
                #(x,y) = numpy.unravel_index(places[index],lhi_center.shape)
                #steps.append((x+self.center_r-center_size,y+self.center_c-center_size))

        self.analyse(steps,ns=__main__.__dict__.get('number_sizes',10))
        
            
    def analyse(self,cords,ns=10,absolute=True):
        for (xindex,yindex) in cords:
		if absolute==False:
			xindex = self.center_c + xindex  
			yindex = self.center_r + yindex

                xcoor,ycoor = self.sheet.matrixidx2sheet(xindex,yindex)
                print "Starting surround analysis for cell with index coords and sheet coords: [%d,%d] [%f,%f]"  % (xindex,yindex,xcoor,ycoor) 
                
                c= topo.command.pylabplot.measure_size_response.instance(sheet=self.sheet,num_phase=__main__.__dict__.get('NUM_PHASE',8),num_sizes=ns,max_size=__main__.__dict__.get('MAX_SIZE',1.2),coords=[(xcoor,ycoor)],duration=__main__.__dict__.get('duration',4.0))
                c.duraton=__main__.__dict__.get('duration', 4.0)
                c(coords=[(xcoor,ycoor)],frequencies=[__main__.__dict__.get('FREQ',2.4)])        
                
                self.data_dict[(xindex,yindex)] = {}
                self.data_dict[(xindex,yindex)]["ST"] = self.calculate_RF_sizes(xindex, yindex)
                self.plot_size_tunning(xindex,yindex)
                import pylab
                #pylab.show()
                self.data_dict[(xindex,yindex)]["OCT"] = self.perform_orientation_contrast_analysis(self.data_dict[(xindex,yindex)]["ST"],xcoor,ycoor,xindex,yindex)
                self.plot_orientation_contrast_tuning(xindex,yindex)
                
        f = open(normalize_path("dict.dat"),'wb')
        import pickle
        pickle.dump((self.OR,self.OS,self.MR,self.data_dict),f)
        f.close()

        self.plot_average_oct(independent=True)
        self.plot_average_size_tuning_curve(independent=True)
        
        if False:
            self.lhi = compute_local_homogeneity_index(self.sheet.sheet_views['OrientationPreference'].view()[0]*pi,2.0)                
            f = open(normalize_path('lhi2.0.pickle'),'wb')            
            pickle.dump(self.lhi,f)
            f.close()
        else:        
            f = open(normalize_path('lhi2.0.pickle'),'rb')            
            self.lhi = pickle.load(f)

        # determine pinwheels and domain centers
        pinwheels = []
        centers = []
        for coords in self.data_dict.keys():    
            if self.lhi[coords] < 0.5:
               pinwheels.append(coords) 
            if self.lhi[coords] > 0.5:
               centers.append(coords) 
               
        self.plot_average_oct(keys=pinwheels,independent=True,string="pinwheels")
        self.plot_average_oct(keys=centers,independent=True,string="domains")

        
        pylab.figure()
        pylab.imshow(self.lhi)
        pylab.colorbar()
        release_fig("LHI")
        raster_plots_lc,raster_plots_hc = self.plot_map_feature_to_surround_modulation_feature_correlations(self.lhi,"Local Homogeneity Index")

	csis=[v["ORTC"]["info"]["csi"] for v in self.data_dict.values]
	
	pylab.figure()
	pylab.hist(csis)
	release_fig(filename="CSI")
		
        self.orrelations_figure(raster_plots_lc)        
        

    def perform_orientation_contrast_analysis(self,data,xcoor,ycoor,xindex,yindex):

        if __main__.__dict__.get('ContrastCenter','HC') == 'LC':	
                contrast_center = self.low_contrast
                curve = data["Contrast = " + str(self.low_contrast) + "%" ]		
        else:
                contrast_center = self.high_contrast
                curve = data["Contrast = " + str(self.high_contrast) + "%" ]

        curve_data={}

        topo.command.pylabplot.measure_or_tuning(num_phase=__main__.__dict__.get('NUM_PHASE',8),num_orientation=__main__.__dict__.get('num_orientations',12),size=curve["measures"]["peak_near_facilitation"]+__main__.__dict__.get('INC',0.0),curve_parameters=[{"contrast":contrast_center}],display=True,coords=[(xcoor,ycoor)],frequencies=[__main__.__dict__.get('FREQ',2.4)],duration=__main__.__dict__.get('duration',4.0))
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",coords=[(xcoor,ycoor)])

        curve_name_ort = "Contrast = " + str(contrast_center) + "%";

        ar = []
        ors = []


        for o in self.sheet.curve_dict['orientation'][curve_name_ort].keys():
            ar.append(self.sheet.curve_dict['orientation'][curve_name_ort][o].view()[0][xindex][yindex])
            ors.append(o)
        
	if __main__.__dict__.get('Max',False):    
        	peak_or_response = max(ar)
		orr = ors[numpy.argmax(ar)]	
	else:
		#lets find the weighted average orientation preference
		o = numpy.angle(numpy.mean(numpy.multiply(numpy.exp(1j*numpy.array(ors)*2),ar)))/2.0 % numpy.pi
		csi = numpy.abs(numpy.mean(numpy.multiply(numpy.exp(1j*numpy.array(ors)*2),ar)))/numpy.sum(numpy.abs(ar))
		assert (o <= numpy.pi and o >= 0) ,  'wrong angle: %f' % (o)
		orr = ors[numpy.argmin(circular_dist(numpy.array(ors) % numpy.pi,o % numpy.pi,numpy.pi))]
		peak_or_response = ar[ors.index(orr)]

        if __main__.__dict__.get('OrrFullfield',False):
            orr = self.sheet.sheet_views['OrientationPreference'].view()[0][xindex][yindex]*pi



        curve_data["ORTC"]={}
        curve_data["ORTC"]["info"]={}
        curve_data["ORTC"]["info"]["pref_or"]=orr
	
	if not __main__.__dict__.get('Max',False):
		curve_data["ORTC"]["info"]["csi"]=csi

        print "ORIENTATION:", orr 
        
        if __main__.__dict__.get('LimitedSurround',True):
            topo.command.pylabplot.measure_orientation_contrast(sizecenter=curve["measures"]["peak_near_facilitation"],
                                                             orientation_center=orr,
                                                             #phasecenter=self.sheet.sheet_views['PhasePreference'].view()[0][xindex][yindex]*2*pi,
                                                             sizesurround=(curve["measures"]["peak_supression"]+curve["measures"]["peak_near_facilitation"])/2,
                                                             size=0.0,
                                                             display=True,
                                                             contrastcenter=contrast_center,
                                                             thickness=(curve["measures"]["peak_supression"]-curve["measures"]["peak_near_facilitation"])/2,
                                                             duration=__main__.__dict__.get('duration',4.0),
                                                             num_phase=__main__.__dict__.get('NUM_PHASE',8),
                                                             frequencies=[__main__.__dict__.get('FREQ',2.4)],
                                                             curve_parameters=[{"contrastsurround":contrast_center}],coords=[(xcoor,ycoor)])                  
        else:
            topo.command.pylabplot.measure_orientation_contrast(sizecenter=curve["measures"]["peak_near_facilitation"]+__main__.__dict__.get('INC',0.0),
                                                             orientation_center=orr,
                                                             #phasecenter=self.sheet.sheet_views['PhasePreference'].view()[0][xindex][yindex]*2*pi,
                                                             sizesurround=4.0,
                                                             size=0.0,
                                                             display=True,
                                                             contrastcenter=contrast_center,
                                                             thickness=4.0-curve["measures"]["peak_near_facilitation"]-__main__.__dict__.get('SPACE',0.0)-__main__.__dict__.get('INC',0.0),
                                                             duration=__main__.__dict__.get('duration',4.0),
                                                             num_phase=__main__.__dict__.get('NUM_PHASE',8),
				                             frequencies=[__main__.__dict__.get('FREQ',2.4)],
                                                             curve_parameters=[{"contrastsurround":contrast_center}],coords=[(xcoor,ycoor)])
        
	for curve_label in sorted(self.sheet.curve_dict['orientationsurround'].keys()):
   	    # due to rounding errors lets find the orthogonal orientation
	    surr_ors =  numpy.array(self.sheet.curve_dict['orientationsurround'][curve_label].keys())
            orr_ort = surr_ors[numpy.argmin(numpy.abs(surr_ors - orr - numpy.pi/2.0))]
            orr = surr_ors[numpy.argmin(numpy.abs(surr_ors - orr))]
		
            print curve_label
            curve_data[curve_label]={}
            curve_data[curve_label]["data"]=self.sheet.curve_dict['orientationsurround'][curve_label]
            curve_data[curve_label]["measures"]={}
            print self.sheet.curve_dict['orientationsurround'][curve_label].keys() , "\nAAA" , orr_ort," ", orr
            pref_or_resp=self.sheet.curve_dict['orientationsurround'][curve_label][orr].view()[0][xindex][yindex]
            cont_or_resp=self.sheet.curve_dict['orientationsurround'][curve_label][orr_ort].view()[0][xindex][yindex]
            
            
            if peak_or_response != 0:
                curve_data[curve_label]["measures"]["or_suppression"]=(cont_or_resp-pref_or_resp)/peak_or_response
            else: 
                curve_data[curve_label]["measures"]["or_suppression"]=0.0
            
            x_values = sorted(curve_data[curve_label]["data"].keys())
            y_values = []
            for k in x_values:
                y_values.append(curve_data[curve_label]["data"][k].view()[0][xindex, yindex])

            ssi = self.calculate_octc_selectivity(2*x_values,peak_or_response-y_values)
            curve_data[curve_label]["measures"]["SSI"]=ssi



        curve_name_orc = "Contrastsurround = " + str(contrast_center) + "%";

        pref_or_resp=self.sheet.curve_dict['orientationsurround'][curve_name_orc][orr].view()[0][xindex][yindex]
        cont_or_resp=self.sheet.curve_dict['orientationsurround'][curve_name_orc][orr_ort].view()[0][xindex][yindex]


        curve_data["ORTC"]["data"]=self.sheet.curve_dict['orientation'][curve_name_ort]
        curve_data["ORTC"]["measures"]={}
        curve_data["ORTC"]["measures"]["colinear_lc_suppresion_index"] = (peak_or_response - pref_or_resp) / peak_or_response
        curve_data["ORTC"]["measures"]["orcontrast_lc_suppresion_index"] = (peak_or_response - cont_or_resp) / peak_or_response


        return curve_data 


    
    def calculate_RF_sizes(self,xindex, yindex):
        curve_data = {}
        hc_curve_name = "Contrast = " + str(self.high_contrast) + "%";
        lc_curve_name = "Contrast = " + str(self.low_contrast) + "%";
        for curve_label in [hc_curve_name,lc_curve_name]:
            curve = self.sheet.curve_dict['size'][curve_label]
            curve_data[curve_label] = {}
            curve_data[curve_label]["data"] = curve

            x_values = sorted(curve.keys())
            y_values = [curve[key].view()[0][xindex, yindex] for key in x_values]

            #compute critical indexes in the size tuning curves
            curve_data[curve_label]["measures"]={}
            curve_data[curve_label]["measures"]["peak_near_facilitation_index"] = numpy.argmax(y_values)
            curve_data[curve_label]["measures"]["peak_near_facilitation"] = x_values[curve_data[curve_label]["measures"]["peak_near_facilitation_index"]]

            if(curve_data[curve_label]["measures"]["peak_near_facilitation_index"] < (len(y_values) - 1)):
		print y_values
		print curve_label
                curve_data[curve_label]["measures"]["peak_supression_index"] = curve_data[curve_label]["measures"]["peak_near_facilitation_index"] + numpy.argmin(y_values[curve_data[curve_label]["measures"]["peak_near_facilitation_index"] + 1:]) + 1
                curve_data[curve_label]["measures"]["peak_supression"] = x_values[curve_data[curve_label]["measures"]["peak_supression_index"]]
                curve_data[curve_label]["measures"]["suppresion_index"] = (y_values[curve_data[curve_label]["measures"]["peak_near_facilitation_index"]] - y_values[curve_data[curve_label]["measures"]["peak_supression_index"]])/ y_values[curve_data[curve_label]["measures"]["peak_near_facilitation_index"]]
	    else:
	        curve_data[curve_label]["measures"]["peak_supression_index"]=len(y_values)-1
		curve_data[curve_label]["measures"]["peak_supression"]=x_values[curve_data[curve_label]["measures"]["peak_supression_index"]]
		curve_data[curve_label]["measures"]["suppresion_index"]=0

            if(curve_data[curve_label]["measures"].has_key("peak_supression_index") and (curve_data[curve_label]["measures"]["peak_supression_index"] < (len(y_values) - 1))):
                curve_data[curve_label]["measures"]["peak_far_facilitation_index"] = curve_data[curve_label]["measures"]["peak_supression_index"] + numpy.argmax(y_values[curve_data[curve_label]["measures"]["peak_supression_index"] + 1:]) + 1
                curve_data[curve_label]["measures"]["peak_far_facilitation"] = x_values[curve_data[curve_label]["measures"]["peak_far_facilitation_index"]]
                curve_data[curve_label]["measures"]["counter_suppresion_index"] = (y_values[curve_data[curve_label]["measures"]["peak_far_facilitation_index"]] - y_values[curve_data[curve_label]["measures"]["peak_supression_index"]])/ y_values[curve_data[curve_label]["measures"]["peak_near_facilitation_index"]]


        curve_data[hc_curve_name]["measures"]["contrast_dependent_shift"]=curve_data[lc_curve_name]["measures"]["peak_near_facilitation"]/curve_data[hc_curve_name]["measures"]["peak_near_facilitation"]             
        curve_data[lc_curve_name]["measures"]["contrast_dependent_shift"]=curve_data[lc_curve_name]["measures"]["peak_near_facilitation"]/curve_data[hc_curve_name]["measures"]["peak_near_facilitation"]
        return curve_data


    def plot_size_tunning(self, xindex, yindex):
        fig = pylab.figure()
        #f = fig.add_subplot(111, autoscale_on=False, xlim=(-0.1, 3.0), ylim=(-0.1, 4.0))
        f = fig.add_subplot(111)
        pylab.title(self.sheet_name, fontsize=12)
        colors=['red','blue','green','purple','orange','black','yellow']
        
        measurment = self.data_dict[(xindex,yindex)]["ST"]
        i = 0
        for curve_label in measurment.keys():
            curve =  measurment[curve_label]["data"]
            x_values = sorted(curve.keys())
            y_values = [curve[key].view()[0][xindex, yindex] for key in x_values]
            
            f.plot(x_values, y_values, lw=3, color=colors[i],label=curve_label)
            
                
            f.annotate('', xy=(measurment[curve_label]["measures"]["peak_near_facilitation"], y_values[measurment[curve_label]["measures"]["peak_near_facilitation_index"]]), xycoords='data',
            xytext=(-1, 20), textcoords='offset points', arrowprops=dict(facecolor='green', shrink=0.05))
    
    
            if measurment[curve_label]["measures"].has_key("peak_supression"):
                f.annotate('', xy=(measurment[curve_label]["measures"]["peak_supression"], y_values[measurment[curve_label]["measures"]["peak_supression_index"]]), xycoords='data',
                           xytext=(-1, 20), textcoords='offset points', arrowprops=dict(facecolor='red', shrink=0.05))
            
            if measurment[curve_label]["measures"].has_key("peak_far_facilitation"):
                f.annotate('', xy=(measurment[curve_label]["measures"]["peak_far_facilitation"], y_values[measurment[curve_label]["measures"]["peak_far_facilitation_index"]]), xycoords='data',
                           xytext=(-1, 20), textcoords='offset points', arrowprops=dict(facecolor='blue', shrink=0.05))
            i+=1
            
        release_fig("STC[" + str(xindex) + "," + str(yindex) + "]")


    def calculate_octc_selectivity(self,angles,responses):
        c = 0
        n = 0
        for a,r in zip(angles,responses):
            c = c + r * complex(numpy.cos(a),numpy.sin(a))    
            n = n + numpy.abs(r)
        return numpy.abs(c)/n

    def plot_orientation_contrast_tuning(self, xindex, yindex):
        fig = pylab.figure()
        f = fig.add_subplot(111, autoscale_on=True)
        pylab.title(self.sheet_name, fontsize=12)
        colors=['red','blue','green','purple','orange','black','yellow']
        
        orientation = self.data_dict[(xindex,yindex)]["OCT"]["ORTC"]["info"]["pref_or"]
        
        measurment = self.data_dict[(xindex,yindex)]["OCT"]
        i = 0
        for curve_label in measurment.keys():
            curve =  measurment[curve_label]["data"]
            x_values = sorted(curve.keys())
            y_values = []
            for k in x_values:
                y_values.append(curve[k].view()[0][xindex, yindex])
            
            x_values=numpy.array(x_values)-orientation
	    
	    print x_values

            for j in xrange(0,len(x_values)):
                if x_values[j] > numpy.pi/2.0:
                   x_values[j] -= numpy.pi 
                if x_values[j] < -numpy.pi/2.0:
                   x_values[j] += numpy.pi

            for j in xrange(0,len(x_values)):
                if x_values[j] > numpy.pi/2.0:
                   x_values[j] -= numpy.pi 
                if x_values[j] < -numpy.pi/2.0:
                   x_values[j] += numpy.pi


            inds = numpy.argsort(x_values)
	    y_values = numpy.take(y_values, inds)
            x_values = sorted(x_values)

            numpy.append(y_values,y_values[0])
            numpy.append(x_values,x_values[0]+numpy.pi)
            
            f.plot(x_values, y_values, lw=3, color=colors[i],label=curve_label)
            i+=1
        
	pylab.legend(loc='lower left')
        release_fig("OCTC[" + str(xindex) + "," + str(yindex) + "]")
        
        fig = pylab.figure()
        f = fig.add_subplot(111, autoscale_on=True)

        curve =  measurment["ORTC"]["data"]
        x_values= sorted(curve.keys())
        y_values=[curve[key].view()[0][xindex,yindex] for key in x_values]

        f.plot(x_values, y_values, lw=3)
        release_fig("OTC[" + str(xindex) + "," + str(yindex) + "]")

    def plot_average_size_tuning_curve(self,independent=True):
        
        if independent:
            fig = pylab.figure()
                
        average_curves={}        
        sem_curves={}        
        curves={}        
        
        
        
        for k in self.data_dict.keys():
            # find maximum for the curves
            m = []
            for curve_label in self.data_dict[k]["ST"].keys():
                    xindex, yindex = k
                    curve =  self.data_dict[k]["ST"][curve_label]["data"]
                    x_values = sorted(curve.keys())
                    m.append(numpy.max([curve[key].view()[0][xindex, yindex] for key in x_values]))
            
            m = numpy.max(m)
            
            for curve_label in self.data_dict[k]["ST"].keys():
                    xindex, yindex = k
                    curve =  self.data_dict[k]["ST"][curve_label]["data"]
                    x_values = sorted(curve.keys())
                    y_values = [curve[key].view()[0][xindex, yindex] for key in x_values]
                    
                    # normalize curve
                    y_values = y_values / m
      
                    if curves.has_key(curve_label):
                      curves[curve_label].append(numpy.array(y_values))
                    else:
                      curves[curve_label]=[numpy.array(y_values)]
            
        
        
        
        for curve_label in curves:
            average_curves[curve_label]=numpy.mean(numpy.array(curves[curve_label]),axis=0)
            sem_curves[curve_label]=scipy.stats.sem(numpy.array(curves[curve_label]),axis=0)
        
        
        colors=['red','blue']
        i=0
        for curve_label in average_curves:
            pylab.plot(x_values, average_curves[curve_label]*100, lw=3, color=colors[i]) 
            pylab.errorbar(x_values, average_curves[curve_label]*100, lw=1, ecolor=colors[i], yerr=sem_curves[curve_label]*100, fmt=None) 
            pylab.xticks([0,numpy.max(x_values)/2,numpy.max(x_values)])
            pylab.yticks([0,50,100])
            pylab.gca().set_yticklabels(['0%','50%','100%'])
            pylab.xlim(-0.1,numpy.max(x_values)+0.1)
            pylab.ylim(0,100)
            i=i+1

        disable_top_right_axis(pylab.gca())

        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)  
        pylab.gca().set_xticklabels(pylab.gca().get_xticks(), fontProperties)
	pylab.gca().set_yticklabels(pylab.gca().get_yticks(), fontProperties)

        if independent:
            release_fig("AverageSTC")
   
    def plot_average_oct(self,keys=None,independent=True,string=""):
        if independent:
            fig = pylab.figure()
        
        average_curves={}        
        sem_curves={}        
        curves={}    
        
        if keys == None:
           keys = self.data_dict.keys()
            
        for k in keys:
          xindex, yindex = k          
          curve =  self.data_dict[k]["OCT"]["ORTC"]["data"]
          x_values = sorted(curve.keys())
          m = numpy.max([curve[l].view()[0][xindex, yindex] for l in x_values])
            
          for curve_label in self.data_dict[k]["OCT"].keys():
              if curve_label != 'Contrastsurround = 0%':                
                  orientation = self.data_dict[(xindex,yindex)]["OCT"]["ORTC"]["info"]["pref_or"]
                  curve =  self.data_dict[k]["OCT"][curve_label]["data"]
                  x_values = sorted(curve.keys())
                  y_values = [curve[l].view()[0][xindex, yindex] for l in x_values]
                  x_values=numpy.array(x_values)-orientation
                  
                  for j in xrange(0,len(x_values)):
                      if x_values[j] > numpy.pi/2.0:
                        x_values[j] -= numpy.pi 
                      if x_values[j] < -numpy.pi/2.0:
                        x_values[j] += numpy.pi

                  for j in xrange(0,len(x_values)):
                      if x_values[j] > numpy.pi/2.0:
                        x_values[j] -= numpy.pi 
                      if x_values[j] < -numpy.pi/2.0:
                        x_values[j] += numpy.pi


                  inds = numpy.argsort(x_values)
                  y_values = numpy.take(y_values, inds)
                  x_values = numpy.take(x_values, inds)
                  #x_values = sorted(x_values)
                  
                  if (x_values[0]+numpy.pi/2.0) < 0.0001:
                    y_values = numpy.append(y_values,[y_values[0]])
                    x_values = numpy.append(x_values,[numpy.pi/2])
                  else:
                    y_values = numpy.append([y_values[-1]],y_values)
                    x_values = numpy.append([-numpy.pi/2],x_values)
                  
                  # normalize curve
                  y_values = y_values / m
                    
                  if curves.has_key(curve_label):
                      curves[curve_label].append(numpy.array(y_values))
                  else:
                      curves[curve_label]=[numpy.array(y_values)]
        
        for curve_label in curves:
            average_curves[curve_label]=numpy.mean(numpy.array(curves[curve_label]),axis=0)
            sem_curves[curve_label]=scipy.stats.sem(numpy.array(curves[curve_label]),axis=0)


        colors=['red','blue','green','purple','orange','black','yellow']
        i=0
        
        for curve_label in average_curves.keys():
            pylab.plot(x_values, average_curves[curve_label]*100, lw=3, color=colors[i],label=curve_label) 
            ax = pylab.gca()
            ax.errorbar(x_values, average_curves[curve_label]*100, lw=1, ecolor=colors[i], yerr=sem_curves[curve_label]*100, fmt=None) 
            ax.set_xlim(-numpy.pi/2-0.2,numpy.pi/2.0+0.2)  
            ax.set_xticks([-numpy.pi/2,0,numpy.pi/2.0])
            ax.set_yticks([0,50,100])
            ax.set_yticklabels(['0%','50%','100%'])
            ax.set_ylim(0,110)
            
            
            i=i+1
        #pylab.legend(loc='lower left')
        disable_top_right_axis(pylab.gca())
        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)
        pylab.gca().set_xticklabels([r'-$\frac{\pi}{2}$',r'0',r'$\frac{\pi}{2}$'], fontProperties)
	pylab.gca().set_yticklabels(pylab.gca().get_yticks(), fontProperties)

        if independent:
            release_fig("AverageOCTC" + string) 

    def correlations_figure(self,raster_plots):
        
        self.fig = pylab.figure(facecolor='w',figsize=(15, 18),dpi=800)
        gs = gridspec.GridSpec(3,2)
        gs.update(left=0.07, right=0.95, top=0.95, bottom=0.05,wspace=0.2,hspace=0.2)
        
        ax = pylab.subplot(gs[0,0])
        m,b = numpy.polyfit(raster_plots["or_suppression"][1],raster_plots["or_suppression"][0],1)
        correlation,pval = scipy.stats.pearsonr(raster_plots["or_suppression"][1],raster_plots["or_suppression"][0])
        pylab.scatter(raster_plots["or_suppression"][1],raster_plots["or_suppression"][0],s=18, facecolor = 'r',lw = 0)
        pylab.plot(raster_plots["or_suppression"][1],m*numpy.array(raster_plots["or_suppression"][1])+b,'-k',linewidth=2)
        disable_top_right_axis(ax)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        pylab.ylabel('Orientation-contrast suppression', fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)

        
        ax = pylab.subplot(gs[0,1])
        m,b = numpy.polyfit(raster_plots["SSI"][1],raster_plots["SSI"][0],1)
        correlation,pval = scipy.stats.pearsonr(raster_plots["SSI"][1],raster_plots["SSI"][0])
        pylab.scatter(raster_plots["SSI"][1],raster_plots["SSI"][0],s=18, facecolor = 'r',lw = 0)
        pylab.plot(raster_plots["SSI"][1],m*numpy.array(raster_plots["SSI"][1])+b,'-k',linewidth=2)
        disable_top_right_axis(ax)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        pylab.ylabel('Surround selectivity index', fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)

        
        print raster_plots.keys()
        
        ax = pylab.subplot(gs[1,0])
        m,b = numpy.polyfit(raster_plots["suppresion_index"][1],raster_plots["suppresion_index"][0],1)
        correlation,pval = scipy.stats.pearsonr(raster_plots["suppresion_index"][1],raster_plots["suppresion_index"][0])
        pylab.scatter(raster_plots["suppresion_index"][1],raster_plots["suppresion_index"][0],s=18, facecolor = 'r',lw = 0)
        pylab.plot(raster_plots["suppresion_index"][1],m*numpy.array(raster_plots["suppresion_index"][1])+b,'-k',linewidth=2)
        disable_top_right_axis(ax)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        pylab.ylabel('Suppression index', fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)
        

        ax = pylab.subplot(gs[1,1])
        m,b = numpy.polyfit(raster_plots["counter_suppresion_index"][1],raster_plots["counter_suppresion_index"][0],1)
        correlation,pval = scipy.stats.pearsonr(raster_plots["counter_suppresion_index"][1],raster_plots["counter_suppresion_index"][0])
        pylab.scatter(raster_plots["counter_suppresion_index"][1],raster_plots["counter_suppresion_index"][0],s=18, facecolor = 'r',lw = 0)
        pylab.plot(raster_plots["counter_suppresion_index"][1],m*numpy.array(raster_plots["counter_suppresion_index"][1])+b,'-k',linewidth=2)
        disable_top_right_axis(ax)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        pylab.ylabel('Counter suppression index', fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)
        
        
        ax = pylab.subplot(gs[2,0])
        pylab.scatter(self.lhi.ravel(),self.MR.ravel()*2,s=3, facecolor = 'r',lw = 0)
        #xx,z = running_average(self.lhi.ravel(),self.MR.ravel()*2)
        #pylab.plot(xx,z,'k',lw=3.0)
        disable_top_right_axis(ax)
        pylab.xlabel('Local homogeneity index', fontsize=20)
        pylab.ylabel('Modulation ratio', fontsize=20)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)
        
        
        ax = pylab.subplot(gs[2,1])
        pylab.scatter(self.lhi.ravel(),self.OS.ravel(),s=3, facecolor = 'r',lw = 0)
        #xx,z = running_average(self.lhi.ravel(),self.OS.ravel())
        #pylab.plot(xx,z,'k',lw=3.0)           
        disable_top_right_axis(ax)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        pylab.xlabel('Local homogeneity index', fontsize=20)
        pylab.ylabel('Orientation selectivity', fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'xticklabels'), fontsize=20)
        pylab.setp(pylab.getp(pylab.gca(), 'yticklabels'), fontsize=20)

        print scipy.stats.pearsonr(self.lhi.ravel(),self.OS.ravel())        
        print scipy.stats.pearsonr(self.lhi.ravel(),self.MR.ravel()*2)
        
        release_fig("CorrelationFigure") 


    def plot_map_feature_to_surround_modulation_feature_correlations(self,map_feature,map_feature_name):
            
            from numpy import polyfit
            
            raster_plots_lc={}
            raster_plots_hc={}
            for (xcoord,ycoord) in self.data_dict.keys():
                for curve_type in self.data_dict[(xcoord,ycoord)].keys():
                    print curve_type
                    if curve_type == "ST":
                       curve_label = "Contrast"
                    else:
                       curve_label = "Contrastsurround" 
                    
                    
                    if self.data_dict[(xcoord,ycoord)][curve_type].has_key(curve_label + " = " + str(self.high_contrast) + "%"):
                        for measure_name in self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = " + str(self.high_contrast) + "%"]["measures"].keys():
                            print measure_name
                            if not raster_plots_hc.has_key(measure_name):
                                raster_plots_hc[measure_name]=[[],[]]    
                            raster_plots_hc[measure_name][0].append(self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = " + str(self.high_contrast) + "%"]["measures"][measure_name])
                            raster_plots_hc[measure_name][1].append(map_feature[xcoord,ycoord])        
                    
                    if self.data_dict[(xcoord,ycoord)][curve_type].has_key(curve_label + " = " + str(self.low_contrast) + "%"):
                        for measure_name in self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = " + str(self.low_contrast) + "%"]["measures"].keys():
                            print measure_name
                            if not raster_plots_lc.has_key(measure_name):
                                raster_plots_lc[measure_name]=[[],[]]    
                            raster_plots_lc[measure_name][0].append(self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = "  + str(self.low_contrast) + "%"]["measures"][measure_name])
                            raster_plots_lc[measure_name][1].append(map_feature[xcoord,ycoord])        
            
            for key in raster_plots_hc.keys():
                    fig = pylab.figure()
                    f = fig.add_subplot(111)
                    f.set_xlabel(str(key).replace('_',' '))
                    f.set_ylabel(map_feature_name.replace('_',' '))
                    try:
                        #correlation = numpy.corrcoef(raster_plots_hc[key][0],raster_plots_hc[key][1])[0,1]
                        import scipy.stats
                        correlation = scipy.stats.pearsonr(raster_plots_hc[key][0],raster_plots_hc[key][1])[0]
                        pval= scipy.stats.pearsonr(raster_plots_hc[key][0],raster_plots_hc[key][1])[1] 
                    except FloatingPointError:
                          correlation = 0
                          pval = 0
                    m,b = numpy.polyfit(raster_plots_hc[key][0],raster_plots_hc[key][1],1)
                    f.plot(raster_plots_hc[key][0],raster_plots_hc[key][1],'ro')
                    f.plot(raster_plots_hc[key][0],m*numpy.array(raster_plots_hc[key][0])+b,'-k',linewidth=2)
                    release_fig("RasterHC<" + map_feature_name + ","+ key +  " Corr:"+ str(correlation) + '|'+ str(pval) + ">")
                    
            
            for key in raster_plots_lc.keys():
                    fig = pylab.figure()
                    f = fig.add_subplot(111)
                    f.set_xlabel(str(key).replace('_',' '))
                    f.set_ylabel(map_feature_name.replace('_',' '))
                    m,b = numpy.polyfit(raster_plots_lc[key][0],raster_plots_lc[key][1],1)
                    try:
                        #correlation = numpy.corrcoef(raster_plots_lc[key][0],raster_plots_lc[key][1])[0,1]
                        import scipy.stats
                        correlation = scipy.stats.pearsonr(raster_plots_lc[key][0],raster_plots_lc[key][1])[0]
                        pval= scipy.stats.pearsonr(raster_plots_lc[key][0],raster_plots_lc[key][1])[1] 
                    except FloatingPointError:
                          correlation = 0
                          pval = 0
                    f.plot(raster_plots_lc[key][0],raster_plots_lc[key][1],'ro')
                    f.plot(raster_plots_lc[key][0],m*numpy.array(raster_plots_lc[key][0])+b,'-k',linewidth=2)
                    release_fig("RasterLC<" + map_feature_name + ","+ key + " Corr:"+ str(correlation)+ '|'+ str(pval) + ">")
            return (raster_plots_lc, raster_plots_hc)
        
    def plot_histograms_of_measures(self):
        histograms_lc = {} 
        histograms_hc = {}
        for (xcoord,ycoord) in self.data_dict.keys():
            for curve_type in self.data_dict[(xcoord,ycoord)].keys():
                print curve_type
                if curve_type == "ST":
                   curve_label = "Contrast"
                else:
                   curve_label = "Contrastsurround"
                print self.data_dict[(xcoord,ycoord)][curve_type].keys()   
                for measure_name in self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = " + str(self.high_contrast) + "%"]["measures"].keys():
                    if not histograms_hc.has_key(curve_type + "_" + measure_name):
                        histograms_hc[curve_type + "_" + measure_name]=[]
                    histograms_hc[curve_type + "_" + measure_name].append(self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = " + str(self.high_contrast) + "%"]["measures"][measure_name])

                for measure_name in self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = " + str(self.low_contrast) + "%"]["measures"].keys():
                    if not histograms_lc.has_key(curve_type + "_" + measure_name):
                        histograms_lc[curve_type + "_" + measure_name]=[]
                    histograms_lc[curve_type + "_" + measure_name].append(self.data_dict[(xcoord,ycoord)][curve_type][curve_label + " = " + str(self.low_contrast) + "%"]["measures"][measure_name])
                
        for key in histograms_lc.keys():
                if ((len(histograms_lc[key]) != 0) and (len(histograms_hc[key]) != 0)):
                    fig = pylab.figure()
                    pylab.title(self.sheet_name+ " " + "MeanLC: " + str(numpy.mean(histograms_lc[key])) + "+/-" + str(numpy.std(histograms_lc[key])/ (len(histograms_lc[key])*len(histograms_lc[key]))) + "MeanHC: " + str(numpy.mean(histograms_hc[key])) + "+/-" + str(numpy.std(histograms_hc[key])/ (len(histograms_hc[key])*len(histograms_hc[key]))) , fontsize=12)
                    
                    f = fig.add_subplot(111)
                    f.set_xlabel(str(key))
                    f.set_ylabel('#Cells')
                    mmax = numpy.max(numpy.max(histograms_lc[key]),numpy.max(histograms_lc[key]))
                    mmin = numpy.min(numpy.min(histograms_lc[key]),numpy.min(histograms_lc[key]))
                    bins = numpy.arange(mmin-0.01,mmax+0.01,(mmax+0.01-(mmin-0.01))/10.0)
                    f.hist(histograms_lc[key],bins=bins,normed=False)
                    #f.axvline(x=numpy.mean(histograms_lc[key]),linewidth=4, color='r')
                    release_fig("Histogram<" + key + ">")
                    print len(histograms_lc[key])
                    print key + "LC mean :" + str(numpy.mean(histograms_lc[key]))
                    print key + "HC mean :" + str(numpy.mean(histograms_hc[key]))
                else:
                    print "Histogram ", key , " empty!"


def running_average(x,y):
    s = 0.1
    z = []
    xx = []
    for i in numpy.arange(0.01,1.0,0.1):
        w = []
        for (j,h) in zip(x,y):
           if abs(j-i) < 0.05:
              w.append(h)  
        z.append(numpy.mean(w))
        xx.append(i)
    return xx,z
    
def compute_local_homogeneity_index(or_map,sigma):
    (xsize,ysize) = or_map.shape 
    
    lhi = numpy.zeros(or_map.shape) 
    
    for sx in xrange(0,xsize):
        for sy in xrange(0,ysize):
            lhi_current=[0,0]
            for tx in xrange(0,xsize):
                for ty in xrange(0,ysize):
                    lhi_current[0]+=numpy.exp(-((sx-tx)*(sx-tx)+(sy-ty)*(sy-ty))/(2*sigma*sigma))*numpy.cos(2*or_map[tx,ty])
                    lhi_current[1]+=numpy.exp(-((sx-tx)*(sx-tx)+(sy-ty)*(sy-ty))/(2*sigma*sigma))*numpy.sin(2*or_map[tx,ty])
            lhi[sx,sy]= numpy.sqrt(lhi_current[0]*lhi_current[0] + lhi_current[1]*lhi_current[1])/(2*numpy.pi*sigma*sigma)
            
    return lhi

def release_fig(filename=None):
    import pylab        
    pylab.show._needmain=False
    if filename is not None:
       fullname=filename+str(topo.sim.time())+".png"
       pylab.savefig(normalize_path(fullname))
    else:
       pylab.show()



