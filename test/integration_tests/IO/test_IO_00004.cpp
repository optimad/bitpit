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
* Testing basic configuraton parser features.
*/
int subtest_001()
{
    std::cout << "Testing basic configuraton parser features" << std::endl;

    // Reset global parset
    //
    // This is only neededto get customized root name and version number.
    config::reset("bitpit", 1);

    // Read the configuration file
    std::cout << std::endl;
    std::cout << "Read configuration file..." << std::endl;
    config::read("data/configuration.xml");

    // Dump the configuration
    std::cout << std::endl;
    std::cout << "Dump configuration..." << std::endl;
    config::root.dump(std::cout, 1);

    // Access configuration
    std::cout << std::endl;
    std::cout << "Access configuration..." << std::endl;
    std::cout << "  - Section \"first\" has color..." << config::root["first"].get("color") << std::endl;
    std::cout << "  - Section \"second\" has color..." << config::root["second"].get("color") << std::endl;
    std::cout << "  - Section \"first\" has distance..." << config::root.getSection("first").get("distance") << std::endl;
    std::cout << "  - Section \"second\" has y data..." << config::root["second"]["data"].get("y") << std::endl;
    std::cout << "  - Section \"first\" option count..." << config::root.getSection("first").getOptionCount() << std::endl;
    std::cout << "  - Section \"first\" sub-section count..." << config::root.getSection("first").getSectionCount() << std::endl;

    int firstDistanceInt = config::root["first"].get<int>("distance");
    std::cout << "  - Section \"first\" has distance (int)..." << firstDistanceInt << std::endl;

    double firstDistanceDouble = config::root["first"].get<double>("distance");
    std::cout << "  - Section \"first\" has distance (double)..." << firstDistanceDouble << std::endl;

    double secondDataDouble = config::root["second"]["data"].get<double>("y");
    std::cout << "  - Section \"section\" has y data (double)..." << secondDataDouble << std::endl;

    bool firstExistsBool = config::root["first"].get<bool>("exists");
    std::cout << "  - Section \"first\" has distance..." << firstExistsBool << std::endl;

    std::cout << "  - Non existent option with fallback..." << config::root.getSection("first").get("none", "111") << std::endl;
    std::cout << "  - Non existent option with fallback (int)..." << config::root.getSection("first").get<int>("none", 111) << std::endl;
    std::cout << "  - Non existent option with fallback (double)..." << config::root.getSection("first").get<double>("none", 111.111) << std::endl;

    // Modify the configuration
    std::cout << std::endl;
    std::cout << "Modify configuration..." << std::endl;

    config::root["first"].set("color", "purple");
    config::root["second"].set("color", "orange");

    config::root["first"].set("distance", 111);
    config::root["second"].set("distance", 222);

    config::root["first"].set("exists", true);
    config::root["second"].set("exists", true);

    config::root["first"]["data"].set("x", 111.111);
    config::root["second"]["data"].set("y", 222.222);

    // Dump the configuration
    std::cout << std::endl;
    std::cout << "Dump configuration..." << std::endl;
    config::root.dump(std::cout, 1);

    // Write the configuration file
    std::cout << std::endl;
    std::cout << "Write configuration file..." << std::endl;
    config::write("configuration_updated.xml");

    return 0;
}


/*!
* Subtest 002
*
* Testing configuration string parser, to absorb/flush xml contents directly on string.
*/
int subtest_002()
{
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Testing configuration XML stringfication" << std::endl;
    std::cout << std::endl;

    //create a sample config tree.
    std::shared_ptr<bitpit::ConfigStringParser> writeStr = std::make_shared<bitpit::ConfigStringParser>(false, true); //enable xml format, enable multisection.
    
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
    
    std::string bufferXML;
    //write writeStr Content to the string buffer
    writeStr->write(bufferXML);
    
    std::cout<<bufferXML<<std::endl;

    //read buffer in another configstring parser
    std::shared_ptr<bitpit::ConfigStringParser> parseStr = std::make_shared<bitpit::ConfigStringParser>(false, true); //enable xml format, enable multisection.
    parseStr->read(bufferXML);

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

    // Run the subtests
    log::cout() << "Testing configuration parser" << std::endl;

    int status;
    try {
        status = subtest_001();
        status = std::max(status, subtest_002());
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
