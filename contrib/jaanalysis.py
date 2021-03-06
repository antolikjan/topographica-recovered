import numpy 
import __main__


def t_function():
    """
    Basic example of an analysis command for run_batch; users are
    likely to need something similar but highly customized.
    """
    print 'Called complex_analysis_function'
    import topo
    from topo.command.analysis import save_plotgroup
    from topo.base.projection import ProjectionSheet
    from topo.sheet import GeneratorSheet
    from topo.analysis.featureresponses import SinusoidalMeasureResponseCommand,FeatureCurveCommand
    import contrib.jacommands
    exec "from topo.analysis.vision import analyze_complexity" in __main__.__dict__
    
    print 'Analysing'
    
    import matplotlib
    matplotlib.rc('xtick', labelsize=17)
    matplotlib.rc('ytick', labelsize=17)					
			
    print 'Build a list of all sheets worth measuring'
    f = lambda x: hasattr(x,'measure_maps') and x.measure_maps
    measured_sheets = filter(f,topo.sim.objects(ProjectionSheet).values())
    input_sheets = topo.sim.objects(GeneratorSheet).values()
							    
    print 'Set potentially reasonable defaults; not necessarily useful'
    topo.command.analysis.coordinate=(0.0,0.0)
    if input_sheets:    topo.command.analysis.input_sheet_name=input_sheets[0].name
    if measured_sheets: topo.command.analysis.sheet_name=measured_sheets[0].name
    
    FeatureCurveCommand.curve_parameters=[{"contrast":30},{"contrast":50},{"contrast":70},{"contrast":90}]
    
    import numpy
    # reset treshold and desable noise before measuring maps
    #m = numpy.mean(topo.sim["V1Simple"].output_fns[2].t)
    #topo.sim["V1Simple"].output_fns[2].t*=0
    #topo.sim["V1Simple"].output_fns[2].t+=m
    #sc = topo.sim["V1Simple"].output_fns[1].generator.scale
    #topo.sim["V1Simple"].output_fns[1].generator.scale=0.0
    #save_plotgroup("Orientation Preference and Complexity")
    save_plotgroup("Orientation Preference")



def complex_analysis_function():
    """
    Basic example of an analysis command for run_batch; users are
    likely to need something similar but highly customized.
    """
    print 'Called complex_analysis_function'
    import topo
    from topo.command.analysis import save_plotgroup
    from topo.base.projection import ProjectionSheet
    from topo.sheet import GeneratorSheet
    from topo.analysis.featureresponses import SinusoidalMeasureResponseCommand,FeatureCurveCommand
    import contrib.jacommands
    exec "from topo.analysis.vision import analyze_complexity" in __main__.__dict__
    
    print 'Analysing'
    
    import matplotlib
    matplotlib.rc('xtick', labelsize=17)
    matplotlib.rc('ytick', labelsize=17)					
			
    print 'Build a list of all sheets worth measuring'
    f = lambda x: hasattr(x,'measure_maps') and x.measure_maps
    measured_sheets = filter(f,topo.sim.objects(ProjectionSheet).values())
    input_sheets = topo.sim.objects(GeneratorSheet).values()
							    
    print 'Set potentially reasonable defaults; not necessarily useful'
    topo.command.analysis.coordinate=(0.0,0.0)
    if input_sheets:    topo.command.analysis.input_sheet_name=input_sheets[0].name
    if measured_sheets: topo.command.analysis.sheet_name=measured_sheets[0].name
    
    FeatureCurveCommand.curve_parameters=[{"contrast":30},{"contrast":50},{"contrast":70},{"contrast":90}]
    
    import numpy
    # reset treshold and desable noise before measuring maps
    #m = numpy.mean(topo.sim["V1Simple"].output_fns[2].t)
    #topo.sim["V1Simple"].output_fns[2].t*=0
    #topo.sim["V1Simple"].output_fns[2].t+=m
    #sc = topo.sim["V1Simple"].output_fns[1].generator.scale
    #topo.sim["V1Simple"].output_fns[1].generator.scale=0.0
    a = topo.sim["V1Complex"].in_connections[0].strength
    
    SinusoidalMeasureResponseCommand.scale=__main__.__dict__.get("analysis_scale",2.3)
    MeasureResponseCommand.scale=__main__.__dict__.get("analysis_scale",2.3)


    #if((float(topo.sim.time()) >= 5003.0) and (float(topo.sim.time()) < 5004.0)): 
#	topo.sim["V1Complex"].in_connections[0].strength=0
#	SinusoidalMeasureResponseCommand.frequencies=[3.0]

#    if((float(topo.sim.time()) >= 5005.0) and (float(topo.sim.time()) < 5006.0)): 
#	SinusoidalMeasureResponseCommand.frequencies=[3.0]

#    if((float(topo.sim.time()) >= 5006.0) and (float(topo.sim.time()) < 5007.0)): 
#    	topo.sim["V1Complex"].in_connections[0].strength=0
#	SinusoidalMeasureResponseCommand.frequencies=[2.4]

#    if((float(topo.sim.time()) >= 5007.0) and (float(topo.sim.time()) < 5008.0)): 
#	SinusoidalMeasureResponseCommand.frequencies=[2.4]



#    if((float(topo.sim.time()) >= 10002.0) and (float(topo.sim.time()) < 10003.0)): 
#	topo.sim["V1Complex"].in_connections[0].strength=0
#	SinusoidalMeasureResponseCommand.frequencies=[2.4]

#    if((float(topo.sim.time()) >= 10003.0) and (float(topo.sim.time()) < 10004.0)): 
#	topo.sim["V1Complex"].in_connections[0].strength=0
#	SinusoidalMeasureResponseCommand.frequencies=[3.0]

#    if((float(topo.sim.time()) >= 10004.0) and (float(topo.sim.time()) < 10005.0)): 
#	SinusoidalMeasureResponseCommand.frequencies=[2.4]

#    if((float(topo.sim.time()) >= 10005.0) and (float(topo.sim.time()) < 10006.0)): 
#	SinusoidalMeasureResponseCommand.frequencies=[3.0]


    save_plotgroup("Orientation Preference and Complexity")
    save_plotgroup("Activity")

									
    # Plot all projections for all measured_sheets
    for s in measured_sheets:
        for p in s.projections().values():
            save_plotgroup("Projection",projection=p)

    
    if(float(topo.sim.time()) >= 10005.0): 
        print 'Measuring orientations'
        SinusoidalMeasureResponseCommand.frequencies=[2.4]
        topo.command.pylabplot.measure_or_tuning_fullfield.instance(sheet=topo.sim["V1Complex"])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0]",sheet=topo.sim["V1Complex"],coords=[(0,0)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,0.1]",sheet=topo.sim["V1Complex"],coords=[(0.1,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,-0.1]",sheet=topo.sim["V1Complex"],coords=[(0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,0.1]",sheet=topo.sim["V1Complex"],coords=[(-0.1,0.1)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,-0.1]",sheet=topo.sim["V1Complex"],coords=[(-0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.2,0.2]",sheet=topo.sim["V1Complex"],coords=[(0.2,0.2)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.2,-0.2]",sheet=topo.sim["V1Complex"],coords=[(0.2,-0.2)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.2,0.2]",sheet=topo.sim["V1Complex"],coords=[(-0.2,0.2)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.2,-0.2]",sheet=topo.sim["V1Complex"],coords=[(-0.2,-0.2)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0.1]",sheet=topo.sim["V1Complex"],coords=[(0.0,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,-0.1]",sheet=topo.sim["V1Complex"],coords=[(0.0,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,0]",sheet=topo.sim["V1Complex"],coords=[(-0.1,0.0)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,0]",sheet=topo.sim["V1Complex"],coords=[(0.1,-0.0)])()

        topo.command.pylabplot.measure_or_tuning_fullfield.instance(sheet=topo.sim["V1Simple"])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0,0]",sheet=topo.sim["V1Simple"],coords=[(0,0)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0.1,0.1]",sheet=topo.sim["V1Simple"],coords=[(0.1,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0.1,-0.1]",sheet=topo.sim["V1Simple"],coords=[(0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[-0.1,0.1]",sheet=topo.sim["V1Simple"],coords=[(-0.1,0.1)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[-0.1,-0.1]",sheet=topo.sim["V1Simple"],coords=[(-0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0.2,0.2]",sheet=topo.sim["V1Simple"],coords=[(0.2,0.2)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0.2,-0.2]",sheet=topo.sim["V1Simple"],coords=[(0.2,-0.2)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[-0.2,0.2]",sheet=topo.sim["V1Simple"],coords=[(-0.2,0.2)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[-0.2,-0.2]",sheet=topo.sim["V1Simple"],coords=[(-0.2,-0.2)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0,0.1]",sheet=topo.sim["V1Simple"],coords=[(0.0,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0,-0.1]",sheet=topo.sim["V1Simple"],coords=[(0.0,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[-0.1,0]",sheet=topo.sim["V1Simple"],coords=[(-0.1,0.0)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="SimpleORTC[0.1,0]",sheet=topo.sim["V1Simple"],coords=[(0.1,-0.0)])()

    #topo.sim["V1Simple"].output_fns[1].generator.scale=sc
    topo.sim["V1Complex"].in_connections[0].strength = a
    #topo.sim["V1Complex"].in_connections[0].strength=st

def push_pull_analysis_function():
    print 'Push pull complex_analysis_function'
    import topo
    import numpy
    from topo.command.analysis import save_plotgroup
    from topo.base.projection import ProjectionSheet
    from topo.sheet import GeneratorSheet
    from topo.analysis.featureresponses import SinusoidalMeasureResponseCommand,FeatureCurveCommand
    import contrib.jacommands
    from contrib.push_pull.CCLISSOM_push_pull_extra import check_RF_corrleation_vs_connection_weights_correlation
    from param import normalize_path
    exec "from topo.analysis.vision import analyze_complexity" in __main__.__dict__

    SinusoidalMeasureResponseCommand.frequencies=[2.4]
    SinusoidalMeasureResponseCommand.scale=__main__.__dict__.get("analysis_scale",2.3)
    
    print 'Analysing'
    
    import matplotlib
    matplotlib.rc('xtick', labelsize=17)
    matplotlib.rc('ytick', labelsize=17)					
			
    print 'Build a list of all sheets worth measuring'
    f = lambda x: hasattr(x,'measure_maps') and x.measure_maps
    measured_sheets = filter(f,topo.sim.objects(ProjectionSheet).values())
    input_sheets = topo.sim.objects(GeneratorSheet).values()
							    
    print 'Set potentially reasonable defaults; not necessarily useful'
    topo.command.analysis.coordinate=(0.0,0.0)
    if input_sheets:    topo.command.analysis.input_sheet_name=input_sheets[0].name
    if measured_sheets: topo.command.analysis.sheet_name=measured_sheets[0].name
    
    FeatureCurveCommand.curve_parameters=[{"contrast":30},{"contrast":50},{"contrast":70},{"contrast":90}]
    
    save_plotgroup("Orientation Preference and Complexity")
    save_plotgroup("Activity",normalize="Individually")

									
    # Plot all projections for all measured_sheets
    for s in measured_sheets:
        for p in s.projections().values():
            save_plotgroup("Projection",projection=p,density=3.0)

    print 'Starting push pull analysis'	
    #analyse_push_pull_connectivity()
    check_RF_corrleation_vs_connection_weights_correlation()
    print 'Finished push pull analysis'
    return
    if(float(topo.sim.time()) >= 10005.0): 
        print 'Measuring orientations'
        SinusoidalMeasureResponseCommand.frequencies=[2.4]
        topo.command.pylabplot.measure_or_tuning_fullfield.instance(sheet=topo.sim["V1Simple"])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0]",sheet=topo.sim["V1Simple"],coords=[(0,0)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,0.1]",sheet=topo.sim["V1Simple"],coords=[(0.1,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,-0.1]",sheet=topo.sim["V1Simple"],coords=[(0.1,-0.1)])()
        

def complex_surround_analysis_function():

    """
    Basic example of an analysis command for run_batch; users are
    likely to need something similar but highly customized.
    """
    import topo
    from topo.command.analysis import save_plotgroup
    from topo.analysis.featureresponses import SinusoidalMeasureResponseCommand,FeatureCurveCommand
    from topo.base.projection import ProjectionSheet
    from topo.sheet import GeneratorSheet
    from topo.command import save_snapshot
    from param import normalize_path
    import contrib.cc_lesi.connection_analysis
    import contrib.jacommands
	
    import contrib.surround_analysis_new_cleaned
    exec "from topo.analysis.vision import analyze_complexity" in __main__.__dict__

    import matplotlib
    matplotlib.rc('xtick', labelsize=17)
    matplotlib.rc('ytick', labelsize=17)					
					
    SinusoidalMeasureResponseCommand.frequencies=[2.4]    
    SinusoidalMeasureResponseCommand.scale=__main__.__dict__.get("analysis_scale",2.3)
    from topo.analysis.featureresponses import PatternPresenter            
    PatternPresenter.duration=4.0
    import topo.command.pylabplot
    reload(topo.command.pylabplot)

    # Build a list of all sheets worth measuring
    f = lambda x: hasattr(x,'measure_maps') and x.measure_maps
    measured_sheets = filter(f,topo.sim.objects(ProjectionSheet).values())
    input_sheets = topo.sim.objects(GeneratorSheet).values()
    # Set potentially reasonable defaults; not necessarily useful
    topo.command.analysis.coordinate=(0.0,0.0)
    if input_sheets:    topo.command.analysis.input_sheet_name=input_sheets[0].name
    if measured_sheets: topo.command.analysis.sheet_name=measured_sheets[0].name
    save_plotgroup("Orientation Preference and Complexity")
    save_plotgroup("Activity",normalize='Individually')
    # Plot all projections for all measured_sheets
    for s in measured_sheets:
        for p in s.projections().values():
            save_plotgroup("Projection",projection=p)
    
    contrib.cc_lesi.connection_analysis.Analyse_connectivity()  

    if(float(topo.sim.time()) > 6020.0):
        if __main__.__dict__.get("save",False):
                save_snapshot(normalize_path('snapshot.typ'))

 
        #contrib.surround_analysis.run_dynamics_analysis(0.0,0.0,0.7,__main__.__dict__.get("analysis_scale",0.3))
        topo.command.pylabplot.measure_or_tuning_fullfield.instance(sheet=topo.sim["V1Complex"])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0]",sheet=topo.sim["V1Complex"],coords=[(0,0)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,0.1]",sheet=topo.sim["V1Complex"],coords=[(0.1,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,-0.1]",sheet=topo.sim["V1Complex"],coords=[(0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,0.1]",sheet=topo.sim["V1Complex"],coords=[(-0.1,0.1)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,-0.1]",sheet=topo.sim["V1Complex"],coords=[(-0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.2,0.2]",sheet=topo.sim["V1Complex"],coords=[(0.2,0.2)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.2,-0.2]",sheet=topo.sim["V1Complex"],coords=[(0.2,-0.2)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.2,0.2]",sheet=topo.sim["V1Complex"],coords=[(-0.2,0.2)])()    
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.2,-0.2]",sheet=topo.sim["V1Complex"],coords=[(-0.2,-0.2)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0.1]",sheet=topo.sim["V1Complex"],coords=[(0.0,0.1)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,-0.1]",sheet=topo.sim["V1Complex"],coords=[(0.0,-0.1)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,0]",sheet=topo.sim["V1Complex"],coords=[(-0.1,0.0)])()    
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,0]",sheet=topo.sim["V1Complex"],coords=[(0.1,-0.0)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.3,0.3]",sheet=topo.sim["V1Complex"],coords=[(0.3,0.3)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.3,-0.3]",sheet=topo.sim["V1Complex"],coords=[(0.3,-0.3)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.3,0.3]",sheet=topo.sim["V1Complex"],coords=[(-0.3,0.3)])()    
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.3,-0.3]",sheet=topo.sim["V1Complex"],coords=[(-0.3,-0.3)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.24,0.24]",sheet=topo.sim["V1Complex"],coords=[(0.24,0.24)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.24,-0.24]",sheet=topo.sim["V1Complex"],coords=[(0.24,-0.24)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.24,0.24]",sheet=topo.sim["V1Complex"],coords=[(-0.24,0.42)])()    
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.24,-0.24]",sheet=topo.sim["V1Complex"],coords=[(-0.24,-0.24)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0.24]",sheet=topo.sim["V1Complex"],coords=[(0.0,0.24)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,-0.24]",sheet=topo.sim["V1Complex"],coords=[(0.0,-0.42)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.24,0]",sheet=topo.sim["V1Complex"],coords=[(-0.24,0.0)])()    
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.24,0]",sheet=topo.sim["V1Complex"],coords=[(0.24,-0.0)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0.3]",sheet=topo.sim["V1Complex"],coords=[(0.0,0.3)])()
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,-0.3]",sheet=topo.sim["V1Complex"],coords=[(0.0,-0.3)])()
	#topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.3,0]",sheet=topo.sim["V1Complex"],coords=[(-0.3,0.0)])()    
        #topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.3,0]",sheet=topo.sim["V1Complex"],coords=[(0.3,-0.0)])()
	
	#contrib.surround_analysis_new_cleaned.surround_analysis("V1Complex").run_analysis_with_step_grid(4,4,max_curves=__main__.__dict__.get("max_curves",20))
	contrib.surround_analysis_new_cleaned.surround_analysis("V1Complex").run_lhi_informed_analysis(max_curves=__main__.__dict__.get("max_curves",26),center_size=__main__.__dict__.get("center_size",20))
        #contrib.surround_analysis_new.surround_analysis("V1Complex").analyse([(0,0),(3,0),(-3,0),(0,3),(0,-3),(3,3),(3,-3),(-3,3),(-3,-3),(6,0),(-6,0),(0,6),(0,-6),(6,6),(6,-6),(-6,6),(-6,-6)],__main__.__dict__.get("number_sizes",10))
	#contrib.surround_analysis_new.surround_analysis("V1Complex").analyse([(57,57),(53,67),(57,59),(61,63),(53,49),(67,65),(51,67),(67,61),(55,49),(47,59),(63,51)],__main__.__dict__.get("number_sizes",10))

													    
def v2_analysis_function():
    """
    Basic example of an analysis command for run_batch; users are
    likely to need something similar but highly customized.
    """
    import topo
    from topo.command.analysis import save_plotgroup
    from topo.base.projection import ProjectionSheet
    from topo.sheet import GeneratorSheet
    exec "from topo.analysis.vision import analyze_complexity" in __main__.__dict__
    from param import normalize_path

    topo.sim["V1Simple"].measure_maps = True
    topo.sim["V1Complex"].measure_maps = True
    topo.sim["V2"].measure_maps = True
    
    topo.sim["V2"].in_connections[0].strength=4
    
    save_plotgroup("Orientation Preference and Complexity")    

    # Plot all projections for all measured_sheets
    measured_sheets = [s for s in topo.sim.objects(ProjectionSheet).values()
                       if hasattr(s,'measure_maps') and s.measure_maps]
    for s in measured_sheets:
        for p in s.projections().values():
            save_plotgroup("Projection",projection=p)

    save_plotgroup("Activity")
#    topo.sim["V1Simple"].measure_maps = False
#    topo.sim["V1Complex"].measure_maps = False
        
    save_plotgroup("Corner OR Preference")
    from topo.command import save_snapshot
#    save_snapshot(normalize_path('snapshot.typ'))


activity_history=numpy.array([])
def rf_analysis():
    import topo
    import pylab
    import topo.analysis.vision
    import contrib.jacommands
    from topo.command.analysis import save_plotgroup
    from topo.base.projection import ProjectionSheet
    from topo.sheet import GeneratorSheet
    from topo.command.analysis import measure_or_tuning_fullfield, measure_or_pref
    from topo.command.pylabplot import cyclic_tuning_curve
    from param import normalize_path    
    
    if(float(topo.sim.time()) <=20010): 
        save_plotgroup("Orientation Preference")
        save_plotgroup("Activity")
    
        # Plot all projections for all measured_sheets
        measured_sheets = [s for s in topo.sim.objects(ProjectionSheet).values()
                           if hasattr(s,'measure_maps') and s.measure_maps]
        for s in measured_sheets:
            for p in s.projections().values():
                save_plotgroup("Projection",projection=p)

        prefix="WithGC"   
        measure_or_tuning_fullfield()
        s=topo.sim["V1"]
        cyclic_tuning_curve(filename_suffix=prefix,filename="OrientationTC:V1:[0,0]",sheet=s,coords=[(0,0)],x_axis="orientation")
        cyclic_tuning_curve(filename_suffix=prefix,filename="OrientationTC:V1:[0.1,0.1]",sheet=s,coords=[(0.1,0.1)],x_axis="orientation")
        cyclic_tuning_curve(filename_suffix=prefix,filename="OrientationTC:V1:[-0.1,-0.1]",sheet=s,coords=[(-0.1,-0.1)],x_axis="orientation")
        cyclic_tuning_curve(filename_suffix=prefix,filename="OrientationTC:V1:[0.1,-0.1]",sheet=s,coords=[(0.1,-0.1)],x_axis="orientation")
        cyclic_tuning_curve(filename_suffix=prefix,filename="OrientationTC:V1:[-0.1,0.1]",sheet=s,coords=[(-0.1,0.1)],x_axis="orientation")
    else:
        topo.command.activity_history = numpy.concatenate((contrib.jacommands.activity_history,topo.sim["V1"].activity.flatten()),axis=1)    

    if(float(topo.sim.time()) == 20000): 
        topo.sim["V1"].plastic=False
        contrib.jacommands.homeostatic_analysis_function()

    if(float(topo.sim.time()) == 20001): 
        pylab.figure()

def gc_homeo_af():
    import contrib.jsldefs
    import topo.command.pylabplot
    import contrib.jacommands
    from topo.command.analysis import save_plotgroup
    from topo.analysis.featureresponses import FeatureResponses , PatternPresenter, FeatureMaps            
    #FeatureResponses.repetitions=10

    FeatureMaps.selectivity_multiplier=20

    PatternPresenter.duration=0.2
    PatternPresenter.apply_output_fns=False
    import topo.command.pylabplot
    reload(topo.command.pylabplot)

    
    on = topo.sim["LGNOn"].in_connections[0].strength
    off = topo.sim["LGNOff"].in_connections[0].strength
    if __main__.__dict__.get("GC",False):
       topo.sim["LGNOn"].in_connections[0].strength=0
       topo.sim["LGNOff"].in_connections[0].strength=0
    
    contrib.jsldefs.homeostatic_analysis_function()
    topo.command.pylabplot.fftplot(topo.sim["V1"].sheet_views["OrientationPreference"].view()[0],filename="V1ORMAPFFT")
    
    from topo.misc.filepath import normalize_path, application_path    
    from scipy.io import write_array
    import numpy
    write_array(normalize_path(str(topo.sim.time())+"orprefmap.txt"), topo.sim["V1"].sheet_views["OrientationPreference"].view()[0])
    write_array(normalize_path(str(topo.sim.time())+"orselmap.txt"), topo.sim["V1"].sheet_views["OrientationSelectivity"].view()[0])
    topo.sim["LGNOn"].in_connections[0].strength = on
    topo.sim["LGNOff"].in_connections[0].strength = off

    print float(topo.sim.time())
    if(float(topo.sim.time()) > 19002.0): 
	#topo.sim["V1"].output_fns[2].scale=0.0
	save_plotgroup("Position Preference")
	PatternPresenter.duration=1.0
        PatternPresenter.apply_output_fns=True
	import topo.command.pylabplot
        reload(topo.command.pylabplot)
        topo.command.pylabplot.measure_or_tuning_fullfield.instance(sheet=topo.sim["V1"],repetitions=10)(repetitions=10)
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0]",sheet=topo.sim["V1"],coords=[(0,0)])()

        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0]",sheet=topo.sim["V1"],coords=[(0.1,0)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0]",sheet=topo.sim["V1"],coords=[(0.1,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0]",sheet=topo.sim["V1"],coords=[(0,0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,0.1]",sheet=topo.sim["V1"],coords=[(0.1,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,-0.1]",sheet=topo.sim["V1"],coords=[(0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,0.1]",sheet=topo.sim["V1"],coords=[(-0.1,0.1)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,-0.1]",sheet=topo.sim["V1"],coords=[(-0.1,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.2,0.2]",sheet=topo.sim["V1"],coords=[(0.2,0.2)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.2,-0.2]",sheet=topo.sim["V1"],coords=[(0.2,-0.2)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.2,0.2]",sheet=topo.sim["V1"],coords=[(-0.2,0.2)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.2,-0.2]",sheet=topo.sim["V1"],coords=[(-0.2,-0.2)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,0.1]",sheet=topo.sim["V1"],coords=[(0.0,0.1)])()
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0,-0.1]",sheet=topo.sim["V1"],coords=[(0.0,-0.1)])()
	topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[-0.1,0]",sheet=topo.sim["V1"],coords=[(-0.1,0.0)])()    
        topo.command.pylabplot.cyclic_tuning_curve.instance(x_axis="orientation",filename="ORTC[0.1,0]",sheet=topo.sim["V1"],coords=[(0.1,-0.0)])()

    if(float(topo.sim.time()) > 20000.0): 
        topo.sim["V1"].output_fns[1].plastic=False
        contrib.jacommands.measure_histogram(iterations=1000) 	
    

def saver_function():
    from topo.command import save_snapshot
    save_snapshot(normalize_path('snapshot.typ'))

def empty():
    a = 1

def sa():
    import topo
    from topo.command.analysis import save_plotgroup
    from param import normalize_path
    import contrib.jacommands
    import contrib.surround_analysis_new_cleaned
    from topo.analysis.featureresponses import SinusoidalMeasureResponseCommand,FeatureCurveCommand
    reload(contrib.surround_analysis_new_cleaned)
    s =  normalize_path.prefix
    exec "from topo.analysis.vision import analyze_complexity" in __main__.__dict__
    from topo.analysis.featureresponses import FeatureResponses , PatternPresenter, FeatureMaps
    #PatternPresenter.duration=4.0
    normalize_path.prefix = s
    #SinusoidalMeasureResponseCommand.scale=__main__.__dict__.get("analysis_scale",1.0)
    #SinusoidalMeasureResponseCommand.frequencies=[2.4]
    #save_plotgroup("Orientation Preference and Complexity")
    
    contrib.surround_analysis_new_cleaned.surround_analysis("V1Complex").run_lhi_informed_analysis(max_curves=__main__.__dict__.get("max_curves",26),center_size=__main__.__dict__.get("center_size",20),index=__main__.__dict__.get("index",0))
    #contrib.surround_analysis_new_cleaned.surround_analysis("V1Complex").analyse([__main__.__dict__.get("coords",(0,0))],__main__.__dict__.get("number_sizes",15),absolute=False)

