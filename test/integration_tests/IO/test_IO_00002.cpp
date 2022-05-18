/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_IO.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing json parser.
* BEWARE: if bitpit has no JSON support (RapidJSON package not installed)
* this test does nothing and return .
*/
int subtest_001()
{
    std::cout << "Testing json configuraton parser features" << std::endl;

    // Declare a bitpit config parser with multisections enabled.
    //
    config::reset("bitpit", 1, true);

    // Read the configuration file
    std::cout << std::endl;
    std::cout << "Read json configuration file..." << std::endl;

    config::read("data/configuration.json");

    bitpit::GlobalConfigParser & root = config::root;
    // Dump the configuration
    std::cout << std::endl;
    std::cout << "Dump configuration..." << std::endl;
    root.dump(std::cout, 1);

    // Access configuration
    std::cout << std::endl;
    std::cout << "Access configuration..." << std::endl;
    std::cout << "  - Section \"first\" has color..." << root["first"].get("color") << std::endl;
    std::cout << "  - Section \"second\" has color..." << root["second"].get("color") << std::endl;
    std::cout << "  - Section \"first\" has distance..." << root.getSection("first").get("distance") << std::endl;
    std::cout << "  - Section \"second\" has y data..." << root["second"]["data"].get("y") << std::endl;
    std::cout << "  - Section \"first\" option count..." << root.getSection("first").getOptionCount() << std::endl;
    std::cout << "  - Section \"first\" sub-section count..." << root.getSection("first").getSectionCount() << std::endl;

    int firstDistanceInt = root["first"].get<int>("distance");
    std::cout << "  - Section \"first\" has distance (int)..." << firstDistanceInt << std::endl;

    double firstDistanceDouble = root["first"].get<double>("distance");
    std::cout << "  - Section \"first\" has distance (double)..." << firstDistanceDouble << std::endl;

    double secondDataDouble = root["second"]["data"].get<double>("y");
    std::cout << "  - Section \"second\" has y data (double)..." << secondDataDouble << std::endl;

    bool firstExistsBool = root["first"].get<bool>("exists");
    std::cout << "  - Section \"first\" has exists..." << firstExistsBool << std::endl;

    std::cout << "  - Empty dummy option reading ... -" << root.getSection("first").get("dummy")<<" - "<< std::endl;

    std::cout << "  - Non existent option with fallback..." << root.getSection("first").get("none", "111") << std::endl;
    std::cout << "  - Non existent option with fallback (int)..." << root.getSection("first").get<int>("none", 111) << std::endl;
    std::cout << "  - Non existent option with fallback (double)..." << root.getSection("first").get<double>("none", 111.111) << std::endl;


    Config::MultiSection multisections = root.getSections("ObjArray");
    for(auto sec : multisections){
        std::cout<<"  - ObjArray element has address "<< sec->get("Address") <<std::endl;
    }

    // Write the configuration file in XML and json
    std::cout << std::endl;
    std::cout << "Write configuration file..." << std::endl;
    root.write("config_test00007_updated.xml");
    root.write("config_test00007_updated.json");

    return 0;
}

/*!
* Subtest 002
*
* Testing configuration string parser, to absorb/flush json contents directly on string.
*/
int subtest_002()
{
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Testing configuration JSON stringfication" << std::endl;
    std::cout << std::endl;

    //create a sample config tree.
    std::shared_ptr<bitpit::ConfigStringParser> writeStr = std::make_shared<bitpit::ConfigStringParser>(true, true); //enable json format, enable multisection.
    
    //fill the writeStr tree with some data 
    writeStr->set("ActivateOption", "1");
    writeStr->set("Origin", "Albuquerque");
    writeStr->set("Span", "180.0 48.0 46.0"); 
    writeStr->set("Dimension", "5 7 3");
    {
        auto &subsect = writeStr->addSection("UserData");
        subsect.set("Name","Walter");
        subsect.set("Surname","White");
        subsect.set("Nick","Heisenberg");
        subsect.set("PassPhrase", "Say my name!");
    } 
    {
        auto &subsect = writeStr->addSection("UserData");
        subsect.set("Name","Jesse");
        subsect.set("Surname","Pinkman");
        subsect.set("Nick","Captain Cook");
        subsect.set("PassPhrase", "Yo!");
    } 
    
    std::string bufferJSON;
    //write writeStr Content to the string buffer
    writeStr->write(bufferJSON);
    
    std::cout<<bufferJSON<<std::endl;

    //read buffer in another configstring parser
    std::shared_ptr<bitpit::ConfigStringParser> parseStr=std::make_shared<bitpit::ConfigStringParser>(true, true); //enable json format, enable multisection.
    parseStr->read(bufferJSON);

    bool check = parseStr->hasSection("UserData");
    check = check && parseStr->hasOption("Origin");
    check = check && parseStr->getSections("UserData").size() == 2;


    return int(!check);
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    // Initialize the logger
    log::manager().initialize(log::MODE_COMBINE);

#if HAS_RAPIDJSON_LIB
    // Run the subtests
    log::cout() << "Testing json configuration parser" << std::endl;

    int status;
    try {
        status = subtest_001();
        status = std::max(status, subtest_002());
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what()<<std::endl;
        exit(1);
    }
#else
    log::cout() << "Skipping test: configuration parser was compiled without JSON support." << std::endl;
#endif

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
