//include root
#include <TH1F.h>
#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TMath.h>

//include FWLite
#include <DataFormats/FWLite/interface/Event.h>
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/FWLite/interface/AutoLibraryLoader.h>

//include pat::Muon
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <PhysicsTools/FWLite/interface/TFileService.h>

using namespace std;

int main(int argc, char *argv[])
{
    // load framework libraries
    gSystem->Load("libFWCoreFWLite");

    //enable automatic library loading
    AutoLibraryLoader::enable();

    //list of input file names (must be in EDM format)
    std::vector<std::string> inputFiles;
    inputFiles.push_back("input.root");

    //output file name (plain ROOT)
    std::string outputFile("output.root");

    //book an output file
    fwlite::TFileService fs = fwlite::TFileService(outputFile.c_str());

    //create a directory for TTrees in the output file
    TFileDirectory dir = fs.mkdir("trees");
    
    //create an output TTree in trees/Events
    TTree *out_tree = 0;
    out_tree = dir.make<TTree>("Events", "Events");

    //index of event
    int ievt = 0;

    //number of bytes read from file
    long bytes_read = 0;

    //keep time
    TStopwatch *stopwatch = new TStopwatch();
    stopwatch->Start();

    //vector of muons will be put into this Handle
    edm::Handle<std::vector<pat::Muon>> muons;
    //name of the muon collection to get
    edm::InputTag muonSource("selectedPatMuons");
    
    //create an output branch for muons
    float mu_pt = 0.0;
    out_tree->Branch("muon_pt", &mu_pt);

    //file loop
    for (unsigned int iFile = 0; iFile < inputFiles.size(); ++iFile)
    {
        // open input file (can be located on castor)
        std::cout << "Opening file " << inputFiles[iFile] << std::endl;
        TFile *in_file = TFile::Open(inputFiles[iFile].c_str());
        if ( in_file )
        {
            std::cout << "File opened successfully" << std::endl;
            double file_time = stopwatch->RealTime();
            stopwatch->Continue();
            
            //load FWLite events
            fwlite::Event ev(in_file);
            std::cout << "Looping over " << ev.size() << " events" << std::endl;
            
            //Event loop
            for (ev.toBegin(); !ev.atEnd(); ++ev, ++ievt)
            {
                std::cout << "event " << ievt << std::endl;

                //get the current event
                edm::EventBase const &event = ev;
                
                //get the muons from the event
                event.getByLabel(muonSource, muons);
                //loop over muons. since Handle is a 'pointer', need to dereference using *
                for(const pat::Muon& muon : *muons) {
                    //get the muon pt
                    std::cout << "muon " << muon.pt() << std::endl;
                    mu_pt = muon.pt();
                }
                out_tree->Fill();
            }
             
            in_file->Close();
                
            //print some statistics
            file_time = stopwatch->RealTime() - file_time;
            stopwatch->Continue();
            bytes_read += in_file->GetBytesRead();
            float mb_read = (float)(in_file->GetBytesRead()) / (1024.0 * 1024.0);
            std::cout << "Closing file " << in_file->GetPath() << " with " << mb_read << " Mb read, "
                    << mb_read / (float)file_time << " Mb/s" << std::endl;
        }
        else
        {
            //file could not be opened, skip
            std::cerr << "Couldn't open an input file: " << inputFiles[iFile] << std::endl;
            continue;
        }
    }

    return 0;
}
