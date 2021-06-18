#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//#include "RHAna.h"

/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=1)
{
  const char *hstfilename, *sumfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events");
  //Declare an instance of our code class
  RHAna m_selec;
  
  if(sample==1){
    //Add one file to chain. This is the input file.
    chain->Add("DYJetsToLLM50_RunII2018_Skim.root");
    //Set names of output files.
    hstfilename = "hst_DY.root";
    sumfilename = "sum_DY.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
  if(sample==2){
    //Add one file to chain. This is the input file.
    chain->Add("TTDilep_RunIIAutumn18_tree_8.root");
    //Set names of output files.
    hstfilename = "hst_TT.root";
    sumfilename = "sum_TT.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
  if(sample==3){
    //Add one file to chain. This is the input file.
    chain->Add("WZ.root");
    //Set names of output files.
    hstfilename = "hst_WZ.root";
    sumfilename = "sum_WZ.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
  if(sample==4){
    //Add one file to chain. This is the input file.
    chain->Add("H-ZZ.root");
    //Set names of output files.
    hstfilename = "hst_H->ZZ.root";
    sumfilename = "sum_H->ZZ.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
  if(sample==5){
    //Add one file to chain. This is the input file.
    chain->Add("VLLM100.root");
    //Set names of output files.
    hstfilename = "hst_VLLM100.root";
    sumfilename = "sum_VLLM100.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
  
  
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);
}

