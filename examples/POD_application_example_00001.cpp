/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

/**
 * \example POD_application_example_00001.cpp
 * 
 * \brief POD application.
 * This application computes the POD basis using the setup specified through 
 * an xml dictionary.
 * <b>To print the usage</b>: ./POD_application_example_00001 --help \n
 */ 

#if BITPIT_ENABLE_MPI
#include <mpi.h>
#endif

#include "pod.hpp"

using namespace bitpit;

//=================================================================================== //

/*!
 * \enum Verbose
 * \brief set verbosity of message returned by the execution
 */
enum class Verbose{
    QUIET=0 /**< no info returned*/,
            NORMAL=1 /**< medium info returned*/,
            FULL=2 /**< full info returned*/
};

//=================================================================================== //
/*!
 * \struct InfoBitpodPP
 * \brief database of essential information absorbed from the custom arguments
 * 
 * The class store from argv string of the main the following information:
 *  - ditcname: (string) name of the target xml dictionary
 *  - vconsole: (enum VERBOSE) type of message verbosity returned by the application on console.
 *  - vlog: (enum VERBOSE) type of message verbosity returned by the application on log file.
 */
struct InfoBitpodPP{

    std::string dictName;       /**< target name of xml dictionary */
    Verbose vconsole;           /**< type of console verbosity*/
    Verbose vlog;               /**< type of log file verbosity */

    /*! Base constructor*/
    InfoBitpodPP(){
        dictName    = "pod.xml";
        vlog        = Verbose::FULL;
        vconsole    = Verbose::NORMAL;
    }
    /*! Destructor */
    ~InfoBitpodPP(){};

    /*! Copy constructor */
    InfoBitpodPP(const InfoBitpodPP & other){
        *this = other;
    }

    /*!Assignement operator */
    InfoBitpodPP & operator=(const InfoBitpodPP & other){
        dictName = other.dictName;
        vlog = other.vlog;
        vconsole = other.vconsole;
        return *this;
    }
};

//=================================================================================== //

/*!
 * Read arguments argv of the main
 * \return InfoBitpodPP structure filled;
 */
InfoBitpodPP readArguments(int argc, char*argv[] ){
    //reading arguments
    if(argc <2) {
        std::cout<<"Error. Not enough arguments found launching the application"<<std::endl;
        std::cout<<"Please run application_example_00001 --help for a brief guide on how to use it"<<std::endl;
        exit(1);
    }

    std::unordered_set<std::string> input;
    for(int i=1; i<argc; ++i){
        std::string temp = argv[i];
        input.insert(bitpit::utils::string::trim(temp));
    }

    if(input.count("--help")>0){
        std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
        std::cout<<""<<std::endl;
        std::cout<<"    Brief application_example_00001 helper, version 1.0"<<std::endl;
        std::cout<<""<<std::endl;
        std::cout<<""<<std::endl;
        std::cout<<"    This is the executable command for running POD instructions from XML Control Dictionaries"<<std::endl;
        std::cout<<"    The command needs mandatorily a XML dictionary to run. It can return execution info on   "<<std::endl;
        std::cout<<"    both console(screen) and external log file. As further debug info, it can plot optional    "<<std::endl;
        std::cout<<"    results of its execution.   "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    The full list of options for running the application are: "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --help                         : print this helper "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --dict=<dictionary name>       : full path to the target xml dictionary. Mandatory. "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --vconsole=<quiet/normal/full> : print info on the execution on console(screen).        "<<std::endl;
        std::cout<<"                                     full is meant for full debug message verbosity, normal for "<<std::endl;
        std::cout<<"                                     a medium verbosity, quiet shut off messaging on console.   "<<std::endl;
        std::cout<<"                                     Default verbosity is medium.                               "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --vlog=<quiet/normal/full>     : print info on the execution external file bitpit.log.   "<<std::endl;
        std::cout<<"                                     full is meant for full debug message verbosity, normal for "<<std::endl;
        std::cout<<"                                     a medium verbosity, quiet shut off messaging on file.      "<<std::endl;
        std::cout<<"                                     Default verbosity is full.                                 "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    For any problem, bug and malfunction please contact bitpit developers.                       "<<std::endl;
        std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
        exit(1);
    }

    std::unordered_map<int, std::string> keymap;
    keymap[0] = "dict=";
    keymap[1] = "vlog=";
    keymap[2] = "vconsole=";

    std::map<int, std::string> final_map;
    /*! Visit input list and search for each key string  in key map. If an input string positively match a key,
     * clean the entry from key, give her the key map marker and store it in the final map.
     */

    for(auto val: input){
        std::size_t pos = std::string::npos;
        int counter=0;
        while (pos == std::string::npos && counter <3){
            pos = val.find(keymap[counter]);
            ++counter;
        }  

        if(pos != std::string::npos) final_map[counter-1] =val.substr(pos+keymap[counter-1].size()); 
    }

    if(final_map.count(0) < 1) {
        std::cout<<"Error. Not valid xml dictionary found launching the application"<<std::endl;
        std::cout<<"Please run application_example_00001 --help for a brief guide on how to use it"<<std::endl;
        exit(1);
    }


    /*! now assign arguments to InfoBitpodPP class */
    InfoBitpodPP result;
    result.dictName = final_map[0];

    if(final_map.count(1)){
        int check = -1 + int(final_map[1]=="quiet") + 2*int(final_map[1]=="normal") + 3*int(final_map[1]=="full");
        if(check == -1) check = 2;
        result.vlog = static_cast<Verbose>(check);
    }

    if(final_map.count(2)){
        int check = -1 + int(final_map[2]=="quiet") + 2*int(final_map[2]=="normal") + 3*int(final_map[2]=="full");
        if(check == -1) check = 1;
        result.vconsole = static_cast<Verbose>(check);
    }

    return result;
}

//=================================================================================== //

/*!
 * Read the POD section of the xml dictionary with following parameters:
 *
 * ## <B>POD</B> ...declaration of object of class POD (only one object allowed)
 * ## <B>Directory</B>: folder where the pod dumping files and results are written
 * ## <B>Name</B>: specific name identifying the pod case
 * ## <B>MemoryMode</B>: memory mode of the pod object
 * ## <B>RunMode</B>: execute mode of the pod object
 * ## <B>WriteMode</B>: write mode of the pod object
 * ## <B>ReconstructionMode</B>: reconstruction mode of the pod object
 * ## <B>ErrorMode</B>: error mode of the pod object
 * ## <B>StaticMesh</B>: condition to set if the mesh of the database snapshots is the same or not
 * ## <B>Modes</B>: target number of retained modes (minimum between desired n modes and energy level)
 * ## <B>Energy</B>: target energy level to set the number of retained modes (minimum between desired n modes and energy level)
 * ## <B>Mean</B>: condition to unset if the mean is not computed
 * ## <B>ErrorThreshold</B>: target error threshold to define the error map bounding box
 * ## <B>ErrorScalarFields</B>: target scalar fields to be used in error bounding box evaluation 
 * ## <B>ErrorVectorFields</B>: target vector fields to be used in error bounding box evaluation  
 * ## <B>/POD</B>
 *
 * \param[in] slotXML   reference to a Section slot of bitpit::Config class.
 * \param[in] podInst   reference to instantiated pod object
 */
void read_podXML(const bitpit::Config::Section & slotXML, POD & podInst){

    if(slotXML.hasOption("Directory")){
        std::string input = slotXML.get("Directory");
        input = bitpit::utils::string::trim(input);
        std::string temp = ".";
        if(!input.empty())  
            podInst.setDirectory(input);
        else               
            podInst.setDirectory(temp);
    }

    if(slotXML.hasOption("Name")){
        std::string input = slotXML.get("Name");
        input = bitpit::utils::string::trim(input);
        std::string temp = "pod";
        if(!input.empty())  
            podInst.setName(input);
        else                
            podInst.setName(temp);
    }

    if(slotXML.hasOption("MemoryMode")){
        std::string input = slotXML.get("MemoryMode");
        input = bitpit::utils::string::trim(input);
        if(input =="MEMORY_LIGHT")
            podInst.setMemoryMode(POD::MemoryMode::MEMORY_NORMAL);
        else
            podInst.setMemoryMode(POD::MemoryMode::MEMORY_NORMAL);
    }

    if(slotXML.hasOption("RunMode")){
        std::string input = slotXML.get("RunMode");
        input = bitpit::utils::string::trim(input);
        if(input =="RESTORE")
            podInst.setRunMode(POD::RunMode::RESTORE);
        else
            podInst.setRunMode(POD::RunMode::COMPUTE);
    }

    if(slotXML.hasOption("WriteMode")){
        std::string input = slotXML.get("WriteMode");
        input = bitpit::utils::string::trim(input);
        if(input =="NONE")
            podInst.setWriteMode(POD::WriteMode::NONE);
        else if(input =="DEBUG")
            podInst.setWriteMode(POD::WriteMode::DEBUG);
        else
            podInst.setWriteMode(POD::WriteMode::DUMP);
    }

    if(slotXML.hasOption("ReconstructionMode")){
        std::string input = slotXML.get("ReconstructionMode");
        input = bitpit::utils::string::trim(input);
        if(input =="MINIMIZATION")
            podInst.setReconstructionMode(POD::ReconstructionMode::MINIMIZATION);
        else
            podInst.setReconstructionMode(POD::ReconstructionMode::PROJECTION);
    }
    
    if(slotXML.hasOption("ErrorMode")){
        std::string input = slotXML.get("ErrorMode");
        input = bitpit::utils::string::trim(input);
        if(input =="SINGLE")
            podInst.setErrorMode(POD::ErrorMode::SINGLE);
        else if(input =="COMBINED")
            podInst.setErrorMode(POD::ErrorMode::COMBINED);
        else
            podInst.setErrorMode(POD::ErrorMode::NONE);
    }    

    if(slotXML.hasOption("StaticMesh")){
        std::string input = slotXML.get("StaticMesh");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        podInst.setStaticMesh(value);
    }
    
    if(slotXML.hasOption("Mean")){
        std::string input = slotXML.get("Mean");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        podInst.setUseMean(value);
    }    

    if(slotXML.hasOption("Modes")){
        std::string input = slotXML.get("Modes");
        input = bitpit::utils::string::trim(input);
        std::size_t temp = std::numeric_limits<std::size_t>::max();
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        podInst.setModeCount(temp);
    }

    if(slotXML.hasOption("Energy")){
        std::string input = slotXML.get("Energy");
        input = bitpit::utils::string::trim(input);
        double temp = 100;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        podInst.setEnergyLevel(temp);
    }
    
   /*! Error bounding box */     
    if(slotXML.hasOption("ErrorThreshold")){
        std::string input = slotXML.get("ErrorThreshold");
        input = bitpit::utils::string::trim(input);
        double temp = 0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        podInst.setErrorThreshold(temp);
    }   

    std::vector<std::string> tempsf;
    std::vector<std::array<std::string,3>> tempvf;    
    
    if(slotXML.hasOption("targetErrorScalarFields")){
        std::string input = slotXML.get("targetErrorScalarFields");
        input = bitpit::utils::string::trim(input);

        if(!input.empty()){
            std::stringstream ss(input);
            ss>>tempsf;
        }       
    } 

    if(slotXML.hasOption("targetErrorVectorFields")){
        std::string input = slotXML.get("targetErrorVectorFields");
        input = bitpit::utils::string::trim(input);

        if(!input.empty()){
            std::stringstream ss(input);
            ss>>tempvf;
        }

    }  
    podInst.setTargetErrorFields(tempsf,tempvf); 

    log::cout().setPriority(bitpit::log::NORMAL);
    log::cout()<< "Finished reading XML dictionary"<<std::endl;
    log::cout().setPriority(bitpit::log::DEBUG);      
  
    /*! Resume pod modes in logger */
    {
        std::vector<std::string> emode(2, "restore");
        emode[1] = "compute";

        std::vector<std::string> mmode(2, "normal");
        emode[1] = "light";

        std::vector<std::string> wmode(3, "dump");
        wmode[1] = "debug";
        wmode[2] = "none";

        std::vector<std::string> rmode(2, "projection");
        wmode[1] = "minimization";
        
        std::vector<std::string> errmode(3, "combined");
        wmode[1] = "single";
        wmode[2] = "none";        

        log::cout()<< "Resume of arguments in "<< podInst.getName() << " : " << std::endl;
        log::cout()<< "execution mode:  "<<emode[static_cast<int>(podInst.getRunMode())]<<std::endl;
        log::cout()<< "memory mode:  "<<mmode[static_cast<int>(podInst.getMemoryMode())]<<std::endl;
        log::cout()<< "write mode:  "<<wmode[static_cast<int>(podInst.getWriteMode())]<<std::endl;
        log::cout()<< "reconstruction mode:  "<<rmode[static_cast<int>(podInst.getReconstructionMode())]<<std::endl;
        log::cout()<< "error mode:  "<<errmode[static_cast<int>(podInst.getErrorMode())]<<std::endl;
        log::cout()<< "number of modes:        "<<podInst.getModeCount()<<std::endl;
        log::cout()<< "energy level:        "<<podInst.getEnergyLevel()<<std::endl;
        log::cout()<< "error threshold:        "<<podInst.getErrorThreshold()<<std::endl;        
        log::cout()<< " "<<std::endl;
        log::cout()<< " "<<std::endl;
    }
}

//=================================================================================== //

/*!
 * Read the JobControls section of the xml dictionary with following parameters:
 *
 * ## <B>/JobControls</B> ...application run controls
 * ## <B>doPODbasis</B>: condition to enable whole POD basis computation
 * ## <B>doLeave1out</B>: condition to enable error map computation through leave-1-out
 * ## <B>doBoundingBox</B>: condition to enable error bounding box computation
 * ## <B>/JobControls</B>
 *
 * \param[in] slotXML   reference to a Section slot of bitpit::Config class.
 * \param[in] podInst   reference to pod object
 */
std::vector<bool> read_jobControlsXML(const bitpit::Config::Section & slotXML, POD & podInst){

    BITPIT_UNUSED(podInst);

    std::vector<bool> controls;
    controls.resize(3,false);
    
    if(slotXML.hasOption("doPODbasis")){
        std::string input = slotXML.get("doPODbasis");
        input = bitpit::utils::string::trim(input);
        bool temp=true;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        controls[0]=temp;
    }
    else 
        controls[0]=true;
    
    if(slotXML.hasOption("doLeave1out")){
        std::string input = slotXML.get("doLeave1out");
        input = bitpit::utils::string::trim(input);
        bool temp=false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        controls[1]=temp;
    }
    
    if(slotXML.hasOption("doBoundingBox")){
        std::string input = slotXML.get("doBoundingBox");
        input = bitpit::utils::string::trim(input);
        bool temp=false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        controls[2]=temp;
    }

    return controls;

}

//=================================================================================== //

/*!
 * read a podXML dictionary, made as :
 * ##xml header string
 * ##dictionary podXML header with version
 * ## <B>Database</B> ...declaration of all run cases to be used as pod database
 * ##  <B>Case</B>
 * ##  <B>Directory</B> ...path of the case<B>/Directory</B>
 * ##  <B>Name</B> ...name of the case<B>/Name</B>
 * ##  <B>NSnapshots</B> ...number of snapshots of the case<B>/NSnapshots</B>
 * ## <B>/Case</B>
 * ## <B>/Database</B>
 *
 * ## <B>POD</B> ...declaration of object of class POD (only one object allowed)
 * ## <B>Directory</B>: folder where the pod dumping files and results are written
 * ## <B>Name</B>: specific name identifying the pod case
 * ## <B>MemoryMode</B>: memory mode of the pod object
 * ## <B>RunMode</B>: execute mode of the pod object
 * ## <B>WriteMode</B>: write mode of the pod object
 * ## <B>ReconstructionMode</B>: reconstruction mode of the pod object
 * ## <B>ErrorMode</B>: error mode of the pod object
 * ## <B>StaticMesh</B>: condition to set if the mesh of the database snapshots is the same or not
 * ## <B>Modes</B>: target number of retained modes (minimum between desired n modes and energy level)
 * ## <B>Energy</B>: target energy level to set the number of retained modes (minimum between desired n modes and energy level)
 * ## <B>Mean</B>: condition to unset if the mean is not computed
 * ## <B>ErrorThreshold</B>: target error threshold to define the error map bounding box
 * ## <B>ErrorScalarFields</B>: target scalar fields to be used in error bounding box evaluation 
 * ## <B>ErrorVectorFields</B>: target vector fields to be used in error bounding box evaluation  
 * ## <B>/POD</B>
 *
 * ## <B>Reconstruct</B> ...declaration of cases to be reconstructed
 * ##  <B>Case</B>
 * ##  <B>Directory</B> ...path of the case<B>/Directory</B>
 * ##  <B>Name</B> ...name of the case<B>/Name</B>
 * ##  <B>NSnapshots</B> ...number of snapshots of the case to be reconstructed<B>/NSnapshots</B>
 * ##  <B>Stride</B> ...snapshots stride<B>Stride</B>
 * ##  <B>FirstSnapshot</B> ...index of first snapshot<B>Stride</B>
 * ## <B>/Case</B>
 * ## <B>/Reconstruct</B>
 * 
 * ## <B>Leave1out</B> ...declaration of cases not used in the leave1out method
 * ##  <B>Case</B>
 * ##  <B>Directory</B> ...path of the case<B>/Directory</B>
 * ##  <B>Name</B> ...name of the case<B>/Name</B>
 * ##  <B>NSnapshots</B> ...number of snapshots not used<B>/NSnapshots</B>
 * ##  <B>Stride</B> ...snapshots stride<B>Stride</B>
 * ##  <B>FirstSnapshot</B> ...index of first snapshot<B>Stride</B>
 * ## <B>/Case</B>
 * ## <B>/leave1out</B>
 * 
 * ## <B>/JobControls</B> ...application run controls
 * ## <B>doPODbasis</B>: condition to enable whole POD basis computation
 * ## <B>doLeave1out</B>: condition to enable error map computation through leave-1-out
 * ## <B>doBoundingBox</B>: condition to enable error bounding box computation
 * ## <B>/JobControls</B>
 *
 * \param[out] declared and instantiated pod object
 */
std::vector<bool> read_Dictionary(POD & podInst) {

    std::vector<bool> exeFlags;
    exeFlags.resize(3,false);
    
    log::cout().setPriority(bitpit::log::NORMAL);
    log::cout()<< "Currently reading XML dictionary"<<std::endl;
    log::cout().setPriority(bitpit::log::DEBUG);

    if(config::root.hasSection("POD")){

        bitpit::Config::Section & podXML = config::root.getSection("POD");

        read_podXML(podXML, podInst);

        log::cout() << "...Instantiated pod: "<< podInst.getName() << std::endl;

    }else{
        log::cout().setPriority(bitpit::log::NORMAL);
        log::cout()<<"No POD section available in the XML dictionary"<<std::endl;
        log::cout().setPriority(bitpit::log::DEBUG);
    }

    /*! Snapshots database */
    if(config::root.hasSection("Database")){

        bitpit::Config::Section & blockXML = config::root.getSection("Database");

        for(auto & sect : blockXML.getSections()){

            std::string name = ".";
            std::string path = ".";
            std::size_t ns, first = 0;
            std::size_t stride = 1;

            name = sect.second->get("Case");
            path = sect.second->get("Directory");
            std::string nss = sect.second->get("NSnapshots");
            std::stringstream ss(nss);
            ss >> ns;
            if(sect.second->hasOption("Stride")){
                std::string sstride = sect.second->get("Stride");
                std::stringstream ssstride(sstride);
                ssstride >> stride;
            }
            if(sect.second->hasOption("FirstSnapshot")){
                std::string sfirst = sect.second->get("FirstSnapshot");
                std::stringstream ssfirst(sfirst);
                ssfirst >> first;
            }
            for (std::size_t i = first; i < ns; i+=stride){
                pod::SnapshotFile snap;
                snap.directory = path;
                snap.name = name + "." + std::to_string(i);
                podInst.addSnapshot(snap);
            }

        }

    }else{
        log::cout().setPriority(bitpit::log::NORMAL);
        log::cout()<<"No Database section available in the XML dictionary"<<std::endl;
        log::cout().setPriority(bitpit::log::DEBUG);
    }

    /*! Reconstruction database */
    if(config::root.hasSection("Reconstruction")){

        bitpit::Config::Section & blockXML = config::root.getSection("Reconstruction");

        for(auto & sect : blockXML.getSections()){

            std::string name = ".";
            std::string path = ".";
            std::size_t ns, first = 0;
            std::size_t stride = 1;

            name = sect.second->get("Case");
            path = sect.second->get("Directory");
            std::string nss = sect.second->get("NSnapshots");
            std::stringstream ss(nss);
            ss >> ns;
            if(sect.second->hasOption("Stride")){
                std::string sstride = sect.second->get("Stride");
                std::stringstream ssstride(sstride);
                ssstride >> stride;
            }
            if(sect.second->hasOption("FirstSnapshot")){
                std::string sfirst = sect.second->get("FirstSnapshot");
                std::stringstream ssfirst(sfirst);
                ssfirst >> first;
            }

            for (std::size_t i = first; i < ns; i+=stride){
                pod::SnapshotFile snap;
                snap.directory = path;
                snap.name = name + "." + std::to_string(i);
                podInst.addReconstructionSnapshot(snap);
            }
        }
    }

    /*! Leave1out database */
    if(config::root.hasSection("Leave1out")){
        bitpit::Config::Section & blockXML = config::root.getSection("Leave1out"); 
        for(auto & sect : blockXML.getSections()){

            std::string name = ".";
            std::string path = ".";
            std::size_t ns, first = 0;
            std::size_t stride = 1;

            name = sect.second->get("Case");
            path = sect.second->get("Directory");
            std::string nss = sect.second->get("NSnapshots");
            std::stringstream ss(nss);
            ss >> ns;
            
            if(sect.second->hasOption("Stride")){
                std::string sstride = sect.second->get("Stride");
                std::stringstream ssstride(sstride);
                ssstride >> stride;
            }
            if(sect.second->hasOption("FirstSnapshot")){
                std::string sfirst = sect.second->get("FirstSnapshot");
                std::stringstream ssfirst(sfirst);
                ssfirst >> first;
            }

            for (std::size_t i = first; i < ns; i+=stride){
                pod::SnapshotFile snap;
                snap.directory = path;
                snap.name = name + "." + std::to_string(i);
                podInst.removeLeave1outSnapshot(snap);
            }
        }
    }

    if(config::root.hasSection("JobControls")){
        bitpit::Config::Section & jobXML = config::root.getSection("JobControls");
        exeFlags=read_jobControlsXML(jobXML, podInst);
    }else      
        exeFlags[0]=true;  
    
    log::cout().setPriority(bitpit::log::NORMAL);
    log::cout()<< "Finished reading XML dictionary"<<std::endl;
    log::cout().setPriority(bitpit::log::DEBUG);
    
    return exeFlags;

}


// =================================================================================== //

/*! Core of xml handler*/
void podcore(const InfoBitpodPP & info) {
    /*! Set the logger and the verbosity of the output messages in execution*/
    std::string log = "bitpit";
    bitpit::log::cout(log);
    switch(int(info.vconsole)){
    case 1 :
        bitpit::log::setConsoleVerbosity(log::cout(), bitpit::log::Verbosity::NORMAL);
        break;
    case 2 :
        bitpit::log::setConsoleVerbosity(log::cout(), bitpit::log::Verbosity::DEBUG);
        break;
    case 0 :
        bitpit::log::setConsoleVerbosity(log::cout(), bitpit::log::Verbosity::QUIET);
        break;
    default:
        break;
    }

    switch(int(info.vlog)){
    case 1 :
        bitpit::log::setFileVerbosity(log::cout(), bitpit::log::Verbosity::NORMAL);
        break;
    case 2 :
        bitpit::log::setFileVerbosity(log::cout(), bitpit::log::Verbosity::DEBUG);
        break;
    case 0 :
        bitpit::log::setFileVerbosity(log::cout(), bitpit::log::Verbosity::QUIET);
        break;
    default:
        break;

    }

    /*! Print resume args info.*/
    log::cout().setPriority(bitpit::log::NORMAL);
    {
        std::vector<std::string> verb(3, "quiet");
        verb[1] = "normal";
        verb[2] = "full";

        log::cout()<< "Resume of arguments:"<<std::endl;
        log::cout()<< "dictionary:         "<<info.dictName<<std::endl;
        log::cout()<< "console verbosity:  "<<verb[static_cast<int>(info.vconsole)]<<std::endl;
        log::cout()<< "log file verbosity: "<<verb[static_cast<int>(info.vlog)]<<std::endl;
        log::cout()<< " "<<std::endl;
        log::cout()<< " "<<std::endl;
    }
    log::cout().setPriority(bitpit::log::DEBUG);

    bitpit::config::reset("podXML", 1);
    bitpit::config::read(info.dictName);

    POD podInst;
    /*! Force mesh type voloctree */
    podInst.setMeshType(POD::MeshType::VOLOCTREE);

    /*! Read dictionary */
    std::vector<bool> exes;
    exes=read_Dictionary(podInst);

    /*! Execute */
    log::cout().setPriority(bitpit::log::NORMAL);
    log::cout()<<"Execution of pod... ";
    log::cout().setPriority(bitpit::log::DEBUG);
    
    /*! Execute*/
    if (exes[1])
        podInst.leave1out();
    
    if (exes[2])
        podInst.evalErrorBoundingBox();
    
    if (exes[0])
        podInst.run();

    log::cout().setPriority(bitpit::log::NORMAL);
    log::cout()<<"...execution of pod done"<<std::endl;

    log::cout().setPriority(bitpit::log::DEBUG);
}

//=================================================================================== //

int main( int argc, char *argv[] ) {

#if BITPIT_ENABLE_MPI
    MPI_Init(&argc, &argv);

    {
#endif
        try{
            /*! Read the arguments and run*/
            InfoBitpodPP info = readArguments(argc, argv);
            podcore(info);
        }
        catch(std::exception & e){
            std::cout<<"application_example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if BITPIT_ENABLE_MPI
    }

    MPI_Finalize();
#endif

    return 0;
}
