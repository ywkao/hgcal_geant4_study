import ROOT as R
import numpy as np
from array import array
import sys

def landau(h,input_file):
	#define canvas
	can_Breitwigner = R.TCanvas("can_bw","can_bw") 
	#define roo var to fit
	adc = R.RooRealVar("","",0.0,100.0)
	adc.setRange("signal",0.0,100.0)
	#define roodata hist
	data = R.RooDataHist("data","data",R.RooArgList(adc),h)
	#define fame to plot
        frame = adc.frame(R.RooFit.Bins(200),R.RooFit.Title("adc"))
	#plot data on fame
	#h1.Draw()
        data.plotOn(frame,R.RooFit.MarkerSize(0.9))#,R.RooFit.DrawOption("hist"))
	#deine papameter to  dfinet the breit wigner shape
	mean = R.RooRealVar("mean","mean",50.,1.,70.)
        width = R.RooRealVar("width","width",10,0.1,50)
	#deine pdf
        #BW = R.RooBreitWigner ("BW","BW",mtop,mean,width)
	lnd = R.RooLandau("landau","landau",adc,mean,width) # sigma fixed
	#model = BW.asTF(R.RooArgList(mtop))
	#define Normalization
        #Norm = R.RooRealVar("Norm","Norm",2000,1,100000000)
        #define model
        # model = R.RooAddPdf("model","Total Model",R.RooArgList(BW),R.RooArgList(Norm))

        #RooFitResult* res;
        res = lnd.fitTo(data, R.RooFit.SumW2Error(R.kTRUE),R.RooFit.Save())#,R.RooFit.Range("signal"))
	#res = hist.Fit(model.GetName(), 'ILS', '')
        #draw fit on frame
        lnd.plotOn(frame, R.RooFit.Name("BW"))
        lnd.paramOn(frame,R.RooFit.Layout(0.55, 0.85, 0.85))

	pad1 = R.TPad('pad1', 'pad1', 0.0, 0.0, 1.0, 0.990683)
        pad1.SetBottomMargin(0.089)
        pad1.SetTicky()
        pad1.SetTickx()
	pad1.Draw()
	pad1.cd()

        frame.Draw()
	can_Breitwigner.Update()
	can_Breitwigner.Print("Plots/"+hist_name+"_landaufit.png")

	raw_input()

def get_hitogram_from_file(input_file,Dir_name,hist_name):
	#read the file to get the hustogrms
        Filename = R.TFile(input_file,"Read")
	#Get the Dir from the file 
        Dir = Filename.GetDirectory(Dir_name)
        print(Dir)	
	#print(Dir.subdir)
	#get histogram
	h1 = Dir.Get(hist_name)
	#print "mean = ", h1.GetMean()
        #print "sigma = ",h1.GetStdDev()
	if("ADC_" in hist_name ):
                h1.GetXaxis().SetTitle("ADC counts")
		h1.GetXaxis().SetRangeUser(0.0,100.0)
		if("ADC_SimhitE"in hist_name):
			h1.GetYaxis().SetRangeUser(0.0,300.0) 
			h1.GetXaxis().SetRangeUser(0.0,350.0)
			h1.GetXaxis().SetTitle("ADC counts") 
			h1.GetYaxis().SetTitle("Simhit Energy")
			h1.GetYaxis().SetTitleOffset(1.45)
		if("SimhitE_ADC"in hist_name):
			h1.GetYaxis().SetRangeUser(0.0,300.0)
	  		h1.GetXaxis().SetRangeUser(0.0,250.0)		
			h1.GetXaxis().SetTitle("Simhit Energy")
                        h1.GetYaxis().SetTitle("ADC counts")
			h1.GetYaxis().SetTitleOffset(1.5)
		h1.GetXaxis().SetLabelSize(0.04)
                h1.GetYaxis().SetLabelSize(0.04)
		h1.SetLineColor(R.kBlue)
		h1.SetLineWidth(1)
		#h1.SetTitle("")
	elif("DigiOccupancy" in hist_name):
		h1.GetXaxis().SetTitle("  X (cm)")
		h1.GetYaxis().SetTitle("  Y (cm)")
		h1.GetXaxis().SetRangeUser(-120,120.0)
		h1.GetYaxis().SetRangeUser(-120,120.0)
                h1.GetXaxis().SetLabelSize(0.04)
                h1.GetYaxis().SetLabelSize(0.04)
		h1.GetYaxis().SetTitleOffset(1.5)
		h1.GetXaxis().SetTitleOffset(1.5)
	#h1.SetLineWidth(2)
	h1.SetStats(R.kTRUE)
	h1.SetDirectory(0)
	#h1.Draw()
	#R.gROOT.cd()
	#return h1
	#raw_input()
	#c1.Update()
	#c1.Print("Plots/"+hist_name+"_HE_scint_D86Geometry_1k_Muons.png")
	#c1.Print("Plots/"+hist_name+"_new_with_xy_slection.png")
	return h1
def get_hist_info(input_file,Dir_name,hist_name):
	Filename = R.TFile(input_file,"Read")
	#Get the Dir from the file 
        Dir = Filename.GetDirectory(Dir_name)
        print(Dir)	
	#get histogram
	h1 = Dir.Get(hist_name)
	st = h1.FindObject("stats")
	print st.Print()
	print "mean = ", h1.GetRMS()
	print "sigma = ",h1.GetStdDev()
	return h1.GetMean()
	
def shower_maximuma(input_file,Dir_name,hist_name):
	c1 = R.TCanvas("c1","c1",0,0,600,600)
	layer_numaber = array( 'd' )
	mean_cell_numaber_involved = array( 'd' )	
	for i in range(1,47):
		layer_numaber.append(i)
		hist_name="cellsnum_perthick_perlayer_120_"+str(i+49)
		mean_cell_numaber_involved.append(get_hist_info(input_file,Dir_name,hist_name))
	gr = R.TGraph(47,layer_numaber,mean_cell_numaber_involved);
	gr.SetLineColor( 2 )
	gr.SetLineWidth( 4 )
	gr.SetMarkerColor( 4 )
	gr.SetMarkerStyle( 21 )
	gr.SetTitle( 'a simple graph' )
	gr.GetXaxis().SetTitle( 'Layer Number' )
	gr.GetYaxis().SetTitle( 'Mean value of number of cells cluster layer' )
	gr.Draw( 'ACP' )
	c1.Update()
	#c1.Print("Plots/"+hist_name+"shower_maximama_wo_cut.png")
	raw_input()	

	print layer_numaber
	print mean_cell_numaber_involved

def create_stand_alone_plots():
        h1 = []
        input_file="/home/mikumar/t3store3/workarea/HGCAL_Validation/Muon_gun/gean_output/geantoutput_merged.root"
        Dir_name="prodEE_DigiSim/"
        hist_name=""
        c1 = R.TCanvas("c1","c1",0,0,600,600)
        R.TGaxis.SetMaxDigits(3)
	#for wafer in ["120","200","300"]:
	#    print wafer
        for layer in range(1,3):#,"2"]:#,"03","04","05","06","07","08","09","08","09",10]:
                #hist_name="ADC_300mum_layer_"+str(layer)
		hist_name="ADC_SimhitE_120mum_layer_"+str(layer)
                print hist_name
                h1.append(get_hitogram_from_file(input_file,Dir_name,hist_name))
                h1[layer-1].Draw()
		c1.Update()
                c1.Print("Plots/"+hist_name+".png")


def create_Landau_fit():

	input_file="geantoutput_10k_0.root"#"DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root"#"DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO_new_D87_nocuton_num_cell.root"
	Dir_name="prodEE_DigiSim/"
	hist_name=""
	h1 = []
	for i in range(1,26):#,"2"]:#,"03","04","05","06","07","08","09","08","09",10]:
		hist_name="ADC_300mum_layer_"+str(i)
		print hist_name
		h1.append(get_hitogram_from_file(input_file,Dir_name,hist_name))
	temphist= h1[0].Clone()
	temphist.Sumw2()
	for i in range(2,25):
		print i
		temphist.Add(h1[i])
	hist_name="ADC_300mum"
	landau(temphist,hist_name)
	
if __name__ == "__main__":
	#create_Landau_fit()	
	create_stand_alone_plots()
