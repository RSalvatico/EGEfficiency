from ROOT import TFile, TCanvas, TH1D, gStyle, TLegend

gStyle.SetOptStat(0)

fIn = [TFile("Efficiency_Standard.root"),TFile("Efficiency_HoE.root")]

h_pT = [fIn[0].Get("EfficiencyCalculator/num_ele_pt_EB"),fIn[0].Get("EfficiencyCalculator/den_ele_pt_EB"),
    fIn[1].Get("EfficiencyCalculator/num_ele_pt_EB"),fIn[1].Get("EfficiencyCalculator/den_ele_pt_EB"),
    fIn[0].Get("EfficiencyCalculator/num_ele_pt_EE"),fIn[0].Get("EfficiencyCalculator/den_ele_pt_EE"),
    fIn[1].Get("EfficiencyCalculator/num_ele_pt_EE"),fIn[1].Get("EfficiencyCalculator/den_ele_pt_EE")]

h_HoE = [fIn[0].Get("EfficiencyCalculator/num_ele_HoE_EB"),fIn[0].Get("EfficiencyCalculator/den_ele_HoE_EB"),
    fIn[1].Get("EfficiencyCalculator/num_ele_HoE_EB"),fIn[1].Get("EfficiencyCalculator/den_ele_HoE_EB"),
    fIn[0].Get("EfficiencyCalculator/num_ele_HoE_EE"),fIn[0].Get("EfficiencyCalculator/den_ele_HoE_EE"),
    fIn[1].Get("EfficiencyCalculator/num_ele_HoE_EE"),fIn[1].Get("EfficiencyCalculator/den_ele_HoE_EE")]

h_eta = [fIn[0].Get("EfficiencyCalculator/num_ele_eta"),fIn[0].Get("EfficiencyCalculator/den_ele_eta"),
    fIn[1].Get("EfficiencyCalculator/num_ele_eta"),fIn[1].Get("EfficiencyCalculator/den_ele_eta")]

h_occupancy = [fIn[0].Get("EfficiencyCalculator/occupancy_phi_eta_all"),fIn[1].Get("EfficiencyCalculator/occupancy_phi_eta_all")]

h_pT[0].Divide(h_pT[0],h_pT[1],1.,1.,"cl=0.683 b(1,1) mode")
h_pT[2].Divide(h_pT[2],h_pT[3],1.,1.,"cl=0.683 b(1,1) mode")
h_pT[4].Divide(h_pT[4],h_pT[5],1.,1.,"cl=0.683 b(1,1) mode")
h_pT[6].Divide(h_pT[6],h_pT[7],1.,1.,"cl=0.683 b(1,1) mode")
#h_HoE[0].Divide(h_HoE[0],h_HoE[1],1.,1.,"cl=0.683 b(1,1) mode")
#h_HoE[2].Divide(h_HoE[2],h_HoE[3],1.,1.,"cl=0.683 b(1,1) mode")
#h_HoE[4].Divide(h_HoE[4],h_HoE[5],1.,1.,"cl=0.683 b(1,1) mode")
#h_HoE[6].Divide(h_HoE[6],h_HoE[7],1.,1.,"cl=0.683 b(1,1) mode")
h_eta[0].Divide(h_eta[0],h_eta[1],1.,1.,"cl=0.683 b(1,1) mode")
h_eta[2].Divide(h_eta[2],h_eta[3],1.,1.,"cl=0.683 b(1,1) mode")

#hRatios = [h_pT[0],h_pT[2],h_pT[4],h_pT[6],h_HoE[0],h_HoE[2],h_HoE[4],h_HoE[6],h_eta[0],h_eta[2]]
hRatios = [h_pT[0],h_pT[2],h_pT[4],h_pT[6],h_eta[0],h_eta[2]]

canvasList = [TCanvas() for i in range(3)]

leg = [TLegend(0.6,0.3,0.8,0.5) for i in range(3)]
for i in range(3):
    leg[i].SetHeader(" ")
    leg[i].SetFillColor(0)
    leg[i].SetBorderSize(0)
    leg[i].SetLineColor(1)
    leg[i].SetLineStyle(1)
    leg[i].SetLineWidth(1)
    leg[i].SetFillStyle(1001)

i_canvas = 0
for a in range(0,len(hRatios)-1):
    hRatios[a].SetMarkerStyle(20)
    hRatios[a].GetYaxis().SetRangeUser(0.,1.1)
    if a % 2 == 0:
        canvasList[i_canvas].cd()
        hRatios[a].Draw("P")
        hRatios[a+1].SetMarkerStyle(26)
        hRatios[a+1].SetMarkerColor(2)
        hRatios[a+1].GetYaxis().SetRangeUser(0.,1.1)
        hRatios[a+1].Draw("SAME,P")
        leg[i_canvas].AddEntry(hRatios[a],"Standard","lep")
        leg[i_canvas].AddEntry(hRatios[a+1],"Tightened HoE","lep")
        leg[i_canvas].Draw("SAME")
        canvasList[i_canvas].SaveAs(hRatios[a].GetName()+".pdf")
        i_canvas += 1

cExtra = [TCanvas() for f in range(4)]
for c in range(2):
    cExtra[c].cd()
    h_occupancy[c].Draw("COLZ")
    if c == 0:
        cExtra[c].SaveAs("occupancy_HLT.pdf")
    else:
        cExtra[c].SaveAs("occupancy_MYHLT.pdf")


h_HoE[0].SetMarkerStyle(20)
h_HoE[4].SetMarkerStyle(20)
h_HoE[2].SetMarkerStyle(26)
h_HoE[6].SetMarkerStyle(26)
h_HoE[2].SetMarkerColor(2)
h_HoE[6].SetMarkerColor(2)
#print("total",h_HoE[0].Integral())
#for i in range(7,16):
#    h_HoE[0].SetBinContent(i,0.)
#print("new",h_HoE[0].Integral())
h_HoE[0].Scale(1/h_HoE[0].Integral())
h_HoE[0].GetYaxis().SetTitle("a. u.")
h_HoE[4].Scale(1/h_HoE[4].Integral())
h_HoE[2].Scale(1/h_HoE[2].Integral())
h_HoE[6].Scale(1/h_HoE[6].Integral())
cExtra[2].cd()
h_HoE[0].Draw("P")
h_HoE[2].Draw("SAME,P")
cExtra[2].SaveAs("HoE_EB.pdf")
cExtra[3].cd()
h_HoE[4].Draw("P")
h_HoE[6].Draw("SAME,P")
cExtra[3].SaveAs("HoE_EE.pdf")


input()




