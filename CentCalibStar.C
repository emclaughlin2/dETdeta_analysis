void CentCalibStar(){
    //read in HIJINGcentcalib.root
    //TFile *fin = new TFile("HIJINGspectrum.root");
    TFile *fin = new TFile("HIJINGspectrum.root");
    //get hnpart
    TH1F *hnpart = (TH1F*)fin->Get("hnpart");
    //find 0-10% centrality bin, 10-20 % ,20-40%, 40-60%, 60-92%
    //loop over xbins
    int totalentries = hnpart->GetEntries();
    int ncumulative = 0;
    int bound5 = 0;
    int bound10 = 0;
    int bound20 = 0;
    int bound40 = 0;
    int bound60 = 0;
    int bound80 = 0;

    //loop from max to min
    for(int i = hnpart->GetNbinsX(); i > 0; i--){
        ncumulative += hnpart->GetBinContent(i);
        if(ncumulative > 0.05*totalentries && hnpart->GetBinLowEdge(i) > bound5){
            
             bound5 = hnpart->GetBinLowEdge(i);
        
        }
	   if(ncumulative > 0.1*totalentries && hnpart->GetBinLowEdge(i) > bound10){
           
            bound10 = hnpart->GetBinLowEdge(i);
            
        }
        if(ncumulative > 0.2*totalentries && hnpart->GetBinLowEdge(i) > bound20){
           
            bound20 = hnpart->GetBinLowEdge(i);
            
        }
        if(ncumulative > 0.4*totalentries && hnpart->GetBinLowEdge(i) > bound40){
           
            bound40 = hnpart->GetBinLowEdge(i);
            
        }
        if(ncumulative > 0.6*totalentries && hnpart->GetBinLowEdge(i) > bound60){
           
            bound60 = hnpart->GetBinLowEdge(i);
            
        }
        if(ncumulative > 0.8*totalentries && hnpart->GetBinLowEdge(i) > bound80){
           
            bound80 = hnpart->GetBinLowEdge(i);
            
        }

   


    }
    std::cout << "0-5% centrality bin: " << bound5 << std::endl;
    std::cout << "5-10% centrality bin: " << bound10 << std::endl;
    std::cout << "10-20% centrality bin: " << bound20 << std::endl;
    std::cout << "20-40% centrality bin: " << bound40 << std::endl;
    std::cout << "40-60% centrality bin: " << bound60 << std::endl;
    std::cout << "60-80% centrality bin: " << bound80 << std::endl;

    //find the average of each centrality bin
    //loop over xbins
    int nentries = 0;
    float sum = 0;
    float average = 0;
    for(int i = 1; i <= hnpart->GetNbinsX(); i++){
        if(hnpart->GetBinLowEdge(i) > bound5){
            nentries += hnpart->GetBinContent(i);
            sum += hnpart->GetBinContent(i)*hnpart->GetBinCenter(i);
        }
    }
    average = sum/nentries;
    std::cout << "0-5% centrality bin average: " << average << std::endl;
    nentries = 0;
    sum = 0;
    average = 0;
    for(int i = 1; i <= hnpart->GetNbinsX(); i++){
        if(hnpart->GetBinLowEdge(i) > bound10 && hnpart->GetBinLowEdge(i) <= bound5){
            nentries += hnpart->GetBinContent(i);
            sum += hnpart->GetBinContent(i)*hnpart->GetBinCenter(i);
        }
    }
    average = sum/nentries;
    std::cout << "5-10% centrality bin average: " << average << std::endl;
    nentries = 0;
    sum = 0;
    average = 0;
    for(int i = 1; i <= hnpart->GetNbinsX(); i++){
        if(hnpart->GetBinLowEdge(i) > bound20 && hnpart->GetBinLowEdge(i) <= bound10){
            nentries += hnpart->GetBinContent(i);
            sum += hnpart->GetBinContent(i)*hnpart->GetBinCenter(i);
        }
    }
    average = sum/nentries;
    std::cout << "10-20% centrality bin average: " << average << std::endl;
    nentries = 0;
    sum = 0;
    average = 0;
    for(int i = 1; i <= hnpart->GetNbinsX(); i++){
        if(hnpart->GetBinLowEdge(i) > bound40 && hnpart->GetBinLowEdge(i) <= bound20){
            nentries += hnpart->GetBinContent(i);
            sum += hnpart->GetBinContent(i)*hnpart->GetBinCenter(i);
        }
    }
    average = sum/nentries;
    std::cout << "20-40% centrality bin average: " << average << std::endl;
    nentries = 0;
    sum = 0;
    average = 0;
    for(int i = 1; i <= hnpart->GetNbinsX(); i++){
        if(hnpart->GetBinLowEdge(i) > bound60 && hnpart->GetBinLowEdge(i) <= bound40){
            nentries += hnpart->GetBinContent(i);
            sum += hnpart->GetBinContent(i)*hnpart->GetBinCenter(i);
        }
    }
    average = sum/nentries;
    std::cout << "40-60% centrality bin average: " << average << std::endl;
    nentries = 0;
    sum = 0;
    average = 0;
    for(int i = 1; i <= hnpart->GetNbinsX(); i++){
        if(hnpart->GetBinLowEdge(i) > bound80 && hnpart->GetBinLowEdge(i) <= bound60){
            nentries += hnpart->GetBinContent(i);
            sum += hnpart->GetBinContent(i)*hnpart->GetBinCenter(i);
        }
    }
    average = sum/nentries;
    std::cout << "60-80% centrality bin average: " << average << std::endl;


        
        
    



}
