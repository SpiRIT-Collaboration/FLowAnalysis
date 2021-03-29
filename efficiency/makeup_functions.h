
TH1 *make_hist_ratio(TH1* hn, TH1* hd, TString name, TString title="")
{
	TH1* hr=nullptr;
	hn->Print();

	if(hn->InheritsFrom("TH2")) hr = (TH2D*)hn->Clone(name);
	else if(hn->InheritsFrom("TH1")) hr = (TH1D*)hn->Clone(name);
	hr->Divide(hd);
	if(!title.IsNull())hr->SetTitle(title);
	return hr;
}

void setup_axis(TAxis *axis,
		double titleSize=0.035, double titleOffset=1.,
		double labelSize=0.035, double labelOffset=0.005
		)
{
	if(axis==nullptr) return 0;
	axis->SetTitleSize(titleSize);
	axis->SetTitleOffset(titleOffset);
	axis->SetLabelSize(labelSize);
	axis->SetLabelOffset(labelOffset);
}

void setup_hist_line(TH1 *h, Color_t color, int style=1, int width=1);
void setup_hist_marker(TH1 *h, Color_t color, int style=20, double size=1.);
void setup_graph_line(TGraph *g, Color_t color, int style=1, int width=1);
void setup_graph_marker(TGraph *g, Color_t color, int style=20, double size=1.);

void setup_line(TObject *o, Color_t color, int style=1, int width=1)
{
	if(o->InheritsFrom("TH1")){
		TH1* hist;
		if(o->InheritsFrom("TH1D"))      { hist = (TH1D*)o; setup_hist_line(hist,color,style,width); }
		else if(o->InheritsFrom("TH1F")) { hist = (TH1F*)o; setup_hist_line(hist,color,style,width); }
		else if(o->InheritsFrom("TH1I")) { hist = (TH1I*)o; setup_hist_line(hist,color,style,width); }
	}

	if(o->InheritsFrom("TGraph")){
		TGraph* graph;
		if(o->InheritsFrom("TGraphAsymmErrors")) { graph = (TGraphAsymmErrors*)o; setup_graph_line(graph,color,style,width); }
		else if(o->InheritsFrom("TGraphErrors")) { graph = (TGraphErrors*)o;      setup_graph_line(graph,color,style,width); }
		else if(o->InheritsFrom("TGraph"))       { graph = (TGraph*)o;            setup_graph_line(graph,color,style,width); }
	}
}
void setup_marker(TObject *o, Color_t color, int style=20, double size=1.)
{
	if(o->InheritsFrom("TH1")){
		TH1* hist;
		if(o->InheritsFrom("TH1D"))      { hist = (TH1D*)o; setup_hist_marker(hist,color,style,size); }
		else if(o->InheritsFrom("TH1F")) { hist = (TH1F*)o; setup_hist_marker(hist,color,style,size); }
		else if(o->InheritsFrom("TH1I")) { hist = (TH1I*)o; setup_hist_marker(hist,color,style,size); }
	}

	if(o->InheritsFrom("TGraph")){
		TGraph* graph;
		if(o->InheritsFrom("TGraphAsymmErrors")) { graph = (TGraphAsymmErrors*)o; setup_graph_marker(graph,color,style,size); }
		else if(o->InheritsFrom("TGraphErrors")) { graph = (TGraphErrors*)o;      setup_graph_marker(graph,color,style,size); }
		else if(o->InheritsFrom("TGraph"))       { graph = (TGraph*)o;            setup_graph_marker(graph,color,style,size); }
	}
}

void setup_hist_line(TH1 *h, Color_t color, int style=1, int width=1)
{
	h->SetLineColor(color);
	h->SetLineStyle(style);
	h->SetLineWidth(width);
}

void setup_hist_marker(TH1 *h, Color_t color, int style=20, double size=1.)
{
	h->SetMarkerColor(color);
	h->SetMarkerStyle(style);
	h->SetMarkerSize(size);
}


void setup_graph_line(TGraph *g, Color_t color, int style=1, int width=1)
{
	g->SetLineColor(color);
	g->SetLineStyle(style);
	g->SetLineWidth(width);
}

void setup_graph_marker(TGraph *g, Color_t color, int style=20, double size=1.)
{
	g->SetMarkerColor(color);
	g->SetMarkerStyle(style);
	g->SetMarkerSize(size);
}


void canvas_partition(TCanvas *C, const Int_t Nx, const Int_t Ny, 
		Float_t lMargin, Float_t rMargin, 
		Float_t bMargin, Float_t tMargin)
{
	if (!C) return;

	// Setup Pad layout:
	Float_t vSpacing = 0.0;
	Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

	Float_t hSpacing = 0.0;
	Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

	Float_t vposd,vposu,vmard,vmaru,vfactor;
	Float_t hposl,hposr,hmarl,hmarr,hfactor;

	for (Int_t i=0;i<Nx;i++) {

		if (i==0) {
			hposl = 0.0;
			hposr = lMargin + hStep;
			hfactor = hposr-hposl;
			hmarl = lMargin / hfactor;
			hmarr = 0.0;
		} else if (i == Nx-1) {
			hposl = hposr + hSpacing;
			hposr = hposl + hStep + rMargin;
			hfactor = hposr-hposl;
			hmarl = 0.0;
			hmarr = rMargin / (hposr-hposl);
		} else {
			hposl = hposr + hSpacing;
			hposr = hposl + hStep;
			hfactor = hposr-hposl;
			hmarl = 0.0;
			hmarr = 0.0;
		}

		for (Int_t j=0;j<Ny;j++) {

			if (j==0) {
				vposd = 0.0;
				vposu = bMargin + vStep;
				vfactor = vposu-vposd;
				vmard = bMargin / vfactor;
				vmaru = 0.0;
			} else if (j == Ny-1) {
				vposd = vposu + vSpacing;
				vposu = vposd + vStep + tMargin;
				vfactor = vposu-vposd;
				vmard = 0.0;
				vmaru = tMargin / (vposu-vposd);
			} else {
				vposd = vposu + vSpacing;
				vposu = vposd + vStep;
				vfactor = vposu-vposd;
				vmard = 0.0;
				vmaru = 0.0;
			}

			C->cd(0);

			TString name = Form("%s_%d_%d",C->GetName(),i,j);
			TPad *pad = (TPad*) gROOT->FindObject(name);
			if (pad) delete pad;
			pad = new TPad(name,"",hposl,vposd,hposr,vposu);
			pad->SetLeftMargin(hmarl);
			pad->SetRightMargin(hmarr);
			pad->SetBottomMargin(vmard);
			pad->SetTopMargin(vmaru);

			pad->SetFrameBorderMode(0);
			pad->SetBorderMode(0);
			pad->SetBorderSize(0);

			pad->Draw();
		}
	}
}
