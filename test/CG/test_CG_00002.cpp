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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include <cmath>
# include <array>
# include <vector>
# include <chrono>
# include <string>
# include <iostream>
#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

using namespace std;
using namespace bitpit;
using namespace bitpit::CGElem;
/*!
 * One by One testing of CGElem namespace methods.
 * 
 */

/*!
 * \return true, test convertBarycentricToFlagSegment passed
 */
bool testConvertBarycentricToFlagSegment(){
    
    bool check(true);
    check &= (convertBarycentricToFlagSegment({{0.5,0.5}}) == 0);
    check &= (convertBarycentricToFlagSegment({{1.0,0.0}}) == 1);
    check &= (convertBarycentricToFlagSegment({{0.0,1.0}}) == 2);
    return check;
}

/*!
 * \return true, test convertBarycentricToFlagTriangle passed
 */
bool testConvertBarycentricToFlagTriangle(){
    
    bool check(true);
    check &= (convertBarycentricToFlagTriangle({{0.1,0.2,0.7}}) == 0);
    check &= (convertBarycentricToFlagTriangle({{1.0,0.0,0.0}}) == 1);
    check &= (convertBarycentricToFlagTriangle({{0.0,0.4,0.6}}) == -2);
    return check;
}

/*!
 * \return true, test convertBarycentricToFlagSimplex passed
 */
bool testConvertBarycentricToFlagSimplex(){
    
    bool check(true);
    check &= (convertBarycentricToFlagSimplex({{0.1,0.2,0.2,0.5}}) == 0);
    check &= (convertBarycentricToFlagSimplex({{1.0,0.0,0.0,0.0}}) == 1);
    check &= (convertBarycentricToFlagSimplex({{0.0,0.4,0.6,0.0}}) == -2);
    return check;
}

/*!
 * \return true, test reconstructPointFromBarycentricSegment passed
 */
bool testReconstructPointFromBarycentricSegment(){
    array3D p1 = {{0.0,0.0,0.0}};
    array3D p2 = {{1.0,1.0,-1.0}};
    array3D test = {{0.75,0.75,-0.75}};
    std::array<double,2> ll = {{0.25,0.75}};
    array3D target = bitpit::CGElem::reconstructPointFromBarycentricSegment(p1, p2,ll.data() );
    return (norm2(target - test)< 1.e-12);
}

/*!
 * \return true, test reconstructPointFromBarycentricTriangle passed
 */
bool testReconstructPointFromBarycentricTriangle(){

    bool check(true);

    array3D p1 = {{0.0,0.0,1.0}};
    array3D p2 = {{1.0,0.0,1.0}};
    array3D p3 = {{0.0,1.0,1.0}};
    array3D test1 = {{0.5,0.5,1.0}};
    array3D test2 = {{0.3,0.5,1.0}};
    array3D target1 = reconstructPointFromBarycentricTriangle(p1, p2, p3, {{0.0,0.5,0.5}});
    array3D target2 = reconstructPointFromBarycentricTriangle(p1, p2, p3, {{0.2,0.3,0.5}});

    check &= (norm2(target1 - test1)< 1.e-12);
    check &= (norm2(target2 - test2)< 1.e-12);
    return check;
}

/*!
 * \return true, test reconstructPointFromBarycentricSimplex passed
 */
bool testReconstructPointFromBarycentricSimplex(){
    std::vector<array3D> pp(5,{{0.0,0.0,0.0}});
    pp[0][2] = 1.0;
    pp[1][0] = 1.0; pp[1][2] = 1.0;
    pp[2][1] = 1.0;
    pp[3][0] = 0.5; pp[3][1] = 0.5;
    std::vector<double> lambda(5,0.2);
    
    array3D test1 = {{0.3,0.3,0.4}};
    array3D target1 = reconstructPointFromBarycentricSimplex(pp, lambda);
    return (norm2(target1 - test1)< 1.e-12);
}

/*!
 * \return true, test projectPointLine passed
 */
bool testProjectPointLine(){

    array3D origin = {{1.0, 2.5, -1.0}};
    array3D dir = {{2.0/7.0, 3.0/7.0, 6.0/7.0}};
    array3D point = {{25.5, 2.5, -1.0}};
    array3D test = {{3.0,5.5,5.0}};
    array3D target = projectPointLine(point, origin, dir);

    return (norm2(target - test)< 1.e-12);
    
}

/*!
 * \return true, test projectPointPlane passed
 */
bool testProjectPointPlane(){
    array3D origin = {{3.0, 6.0, -4.0}};
    array3D normal = {{0,0,1.0}};
    array3D point = {{1.0,1.0,8.0}};
    array3D test = {{1.0,1.0,-4.0}};
    array3D target = projectPointPlane(point, origin, normal);
    
    return (norm2(target - test)< 1.e-12);
}

/*!
 * \return true, test projectPointSegment passed
 */
bool testProjectPointSegment(){
    array3D p1 = {{1.0, 0.0, 2.0}};
    array3D p2 = {{3,0,2.0}};
    array3D point = {{2.0,4.0,2.0}};
    array3D test = {{2,0,2}};
    array3D target = projectPointSegment(point, p1, p2);

    return(norm2(target - test)< 1.e-12);
}

/*!
 * \return true, test projectPointTriangle passed
 */
bool testProjectPointTriangle(){
    array3D p1 = {{1.0, 0.0, 2.0}};
    array3D p2 = {{3,0,2.0}};
    array3D p3 = {{2,2,2.0}};
    array3D point = {{2.0,1.0,0.0}};
    array3D test = {{2,1,2}};
    array3D target = projectPointTriangle(point, p1, p2, p3);

    return (norm2(target - test)< 1.e-12);
}

/*!
 * \return true, test restrictPointTriangle passed
 */
bool testRestrictPointTriangle(){
    array3D p1 = {{1.0,1.0,2.0}};
    array3D p2 = {{2.0,1.0,2.0}};
    array3D p3 = {{1.0,2.0,2.0}};
    array3D test = {{1.0,1.0,2.0}};
    array3D lambda = {{1.5,-0.5,0.0}};
    array3D target = restrictPointTriangle(p1, p2, p3, lambda);
    
    return (norm2(target - test)< 1.e-12);
}


/*!
 * \return true, test projectCloudTriangle passed
 */
bool testProjectCloudTriangle(){

     bool check(true);

     array3D p1 = {{1.0,1.0,2.0}};
     array3D p2 = {{2.0,1.0,2.0}};
     array3D p3 = {{1.0,2.0,2.0}};
     
     std::vector<array3D> cloud;
     std::vector<array3D> test;
     std::vector<array3D> lambda;
     
     cloud.push_back({{1.5,0.0,2.0}});
     cloud.push_back({{1.5,1.5,4.0}});
     
     test.push_back({{1.5,1.0,2.0}});
     test.push_back({{1.5,1.5,2.0}});
     
     std::vector<array3D> projected = projectCloudTriangle(cloud, p1, p2, p3, lambda);
     
     check &= (test.size() == projected.size());
     for(size_t i=0; i<test.size(); ++i){
         check &= (norm2(projected[i] - test[i])< 1.e-12);
     }    
     
     return check;
}


/*!
 * \return true, test projectPointSimplex passed
 */
bool testProjectPointSimplex(){
    
    std::vector<array3D>polygon;
    polygon.push_back({{0.0, 0.0, 2.0}});
    polygon.push_back({{1.0, 0.0, 2.0}});
    polygon.push_back({{1.0, 1.0, 2.0}});
    polygon.push_back({{0.5, 2.0, 2.0}});
    polygon.push_back({{0.0, 1.0, 2.0}});
    
    array3D point = {{0.5,0.73,0.0}};
    array3D test = {{0.5,0.73,2.0}};
    array3D target = projectPointSimplex(point, polygon);
    
    return (norm2(target - test)< 1.e-12);
}

/*!
 * \return true, test projectPointCone passed
 */
bool testProjectPointCone(){
    
    double theta = M_PI/6.0;
    array3D apex = {{0.0,0.0,2.0}};
    array3D axis = {{-1.0,0.0,0.0}};
    array3D dir1 = {{-1.0*cos(theta), sin(theta), 0.0}};
    array3D test = apex  + 2.0*dir1;
    array3D dir2 = {{cos(2*theta), sin(2*theta), 0.0}};
    array3D point = test + 1.8*dir2 ;

    array3D target = projectPointCone(point, apex, axis, theta);
    
    return (norm2(target - test)< 1.e-12);
}

/*!
 * \return true, test distancePointLine passed
 */
bool testDistancePointLine(){
    bool check(true);

    array3D origin = {{1.0, 2.5, -1.0}};
    array3D dir = {{2.0/7.0, 3.0/7.0, 6.0/7.0}};
    array3D point = {{25.5, 2.5,-1.0}};
    array3D test = {{3.0,5.5,5.0}};
    array3D xp;
    double distance = norm2(point - test);
    double target = distancePointLine(point, origin, dir, xp);
    
    check &= (std::abs(target - distance)< 1.e-12);
    check &= (norm2(xp - test) < 1.E-12);
    
    return check;
}

/*!
 * \return true, test distancePointPlane passed
 */
bool testDistancePointPlane(){

    bool check(true);

    array3D origin = {{3.0, 6.0, -4.0}};
    array3D normal = {{0,0,1.0}};
    array3D point = {{1.0,1.0,8.0}};
    array3D test = {{1.0,1.0,-4.0}};
    array3D xp;
    double distance = norm2(point - test);
    double target = distancePointPlane(point, origin, normal,xp);
       
    check &= (std::abs(target - distance)< 1.e-12);
    check &= (norm2(xp - test) < 1.E-12);
    
    return check;
}

/*!
 * \return true, test distancePointSegment passed
 */
bool testDistancePointSegment(){

    bool check(true);

    array3D p1 = {{1.0, 0.0, 2.0}};
    array3D p2 = {{3,0,2.0}};
    array3D point = {{2.0,4.0,2.0}};
    array3D test = {{2,0,2}};
    double distance = norm2(point - test);
    
    std::array<double,2> lambda, lambda_test;
    lambda_test.fill(0.5);
    
    double target0 = distancePointSegment(point, p1, p2);
    double target1 = distancePointSegment(point, p1, p2, lambda);
    
    check &= (std::abs(target0 - distance)< 1.e-12);
    check &= (std::abs(target1 - distance)< 1.e-12);
    check &= (norm2(lambda - lambda_test) < 1.E-12);
    
    return check;
}

/*!
 * \return true, test distancePointTriangle passed
 */
bool testDistancePointTriangle(){
    
    bool check(true);

    array3D p1 = {{1.0, 0.0, 2.0}};
    array3D p2 = {{3,0,2.0}};
    array3D p3 = {{2,2,2.0}};
    array3D point = {{2.0,1.0,0.0}};
    array3D test = {{2,1,2}};
    double distance = norm2(point - test);
    
    std::array<double,3> lambda;
    
    double target0 = distancePointTriangle(point, p1, p2, p3);
    double target1 = distancePointTriangle(point, p1, p2, p3, lambda);
    check &= (std::abs(target0 - distance)< 1.e-12);
    check &= (std::abs(target1 - distance)< 1.e-12);
    
    return check;
}

/*!
 * \return true, test distancePointCone passed
 */
bool testDistancePointCone(){
    double theta = M_PI/6.0;
    array3D apex = {{0.0,0.0,2.0}};
    array3D axis = {{-1.0,0.0,0.0}};
    array3D dir1 = {{-1.0*cos(theta), sin(theta), 0.0}};
    array3D test = apex  + 2.0*dir1;
    array3D dir2 = {{cos(2*theta), sin(2*theta), 0.0}};
    array3D point = test + 1.8*dir2 ;
    
    double distance = norm2(point -test);
    double  target = distancePointCone(point, apex, axis, theta);
    
    return (std::abs(target - distance)< 1.e-12);
}

/*!
 * \return true, test distanceCloudTriangle passed
 */
bool testDistanceCloudTriangle(){
    
    bool check(true);

    array3D p1 = {{1.0,1.0,2.0}};
    array3D p2 = {{2.0,1.0,2.0}};
    array3D p3 = {{1.0,2.0,2.0}};
    std::vector<array3D> cloud, test, lambda;
    cloud.push_back({{1.5,0.0,2.0}});
    cloud.push_back({{1.5,1.5,4.0}});
    
    test.push_back({{1.5,1.0,2.0}});
    test.push_back({{1.5,1.5,2.0}});
    
    std::vector<double> distances(2);
    for(int i=0; i<2; ++i)  distances[i] = norm2(cloud[i] - test[i]);
    
    std::vector<array3D> xps;
    
    auto target0 = distanceCloudTriangle(cloud, p1, p2, p3);
    auto target1 = distanceCloudTriangle(cloud, p1, p2, p3, lambda);
    for(int j=0; j<2; ++j){
        check &= (std::abs(target0[j] - distances[j])< 1.e-12);
        check &= (std::abs(target1[j] - distances[j])< 1.e-12);
    }
    
    return check;
}

/*!
 * \return true, test distancePointSimplex passed
 */
bool testDistancePointSimplex(){
    
    bool check(true);

    std::vector<array3D>polygon;
    polygon.push_back({{0.0, 0.0, 2.0}});
    polygon.push_back({{1.0, 0.0, 2.0}});
    polygon.push_back({{1.0, 1.0, 2.0}});
    polygon.push_back({{0.5, 2.0, 2.0}});
    polygon.push_back({{0.0, 1.0, 2.0}});
    
    array3D point = {{0.5,0.73,0.0}};
    array3D test = {{0.5,0.73,2.0}};
    
    double distance = norm2(point - test);
    std::vector<double> lambda;
    
    auto target0 = distancePointSimplex(point, polygon);
    auto target1 = distancePointSimplex(point, polygon, lambda);
    check &= (std::abs(target0 - distance)< 1.e-12);
    check &= (std::abs(target1 - distance)< 1.e-12);
    
    return check;
}

/*!
 * \return true, test distanceCloudSimplex passed
 */
bool testDistanceCloudSimplex(){
    
    bool check(true);

    std::vector<array3D>polygon;
    polygon.push_back({{0.0, 0.0, 2.0}});
    polygon.push_back({{1.0, 0.0, 2.0}});
    polygon.push_back({{1.0, 1.0, 2.0}});
    polygon.push_back({{0.5, 2.0, 2.0}});
    polygon.push_back({{0.0, 1.0, 2.0}});
    
    std::vector<array3D> cloud, test;
    cloud.push_back({{0.5,0.73,0.0}});
    cloud.push_back({{1.0,2.5,2.0}});
    
    test.push_back({{0.5,0.73,2.0}});
    test.push_back({{0.5,2.0,2.0}});
    
    std::vector<double> distances(2);
    for(int i=0; i<2; ++i)  distances[i] = norm2(cloud[i] - test[i]);
                      
    std::vector<std::vector<double> > lambda; 
    
    auto target0 = distanceCloudSimplex(cloud, polygon);
    auto target1 = distanceCloudSimplex(cloud, polygon, lambda);
    for(int j=0; j<2; ++j){
        check &= (std::abs(target0[j] - distances[j])< 1.e-12);
        check &= (std::abs(target1[j] - distances[j])< 1.e-12);
        
    }
    
    return check;
}


/*!
 * \return true, test distanceLineLine passed
 */
bool testDistanceLineLine(){
    
    bool check(true);
    double cc = pow(2,0.5);
    
    array3D p0 = {{0,0,0}};
    array3D d0 = {{1,0,0}};
    
    array3D p1 = {{2,2,4}};
    array3D d1 = {{1.0/cc,1.0/cc,0.0}};
    
    array3D p2 = {{1.0,0.0,3.0}};
    array3D d2 = {{0,0,-1}};
    
    double dist1 = 4.0;
    double dist2 = 0.0;
    
    auto target1 = distanceLineLine(p0,d0,p1,d1);
    auto target2 = distanceLineLine(p0,d0,p2,d2);
    
    check &= (std::abs(target1-dist1)<1.E-12);
    check &= (std::abs(target2-dist2)<1.E-12);
    
    return check;
}


/*!
 * \return true, test intersectLineLine passed
 */
bool testIntersectLineLine(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{0,0,0}};
    array3D d0 = {{1,0,0}};
    
    array3D p1 = {{-2,2,-3}};
    array3D d1 = {{4.0/cc,-2.0/cc,3.0/cc}};
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectLineLine(p0,d0,p1,d1,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectSegmentSegment passed
 */
bool testIntersectSegmentSegment(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{0,0,0}};
    array3D d0 = {{1,0,0}};
    array3D p1 = p0 + 4.0*d0; 
    
    array3D p2 = {{-2,2,-3}};
    array3D d2 = {{4.0/cc,-2.0/cc,3.0/cc}};
    array3D p3 = p2 + 10.0*d2; 
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectSegmentSegment(p0,p1,p2,p3,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectLinePlane passed
 */
bool testIntersectLinePlane(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{-2,2,-3}};
    array3D d0 = {{4.0/cc,-2.0/cc,3.0/cc}};

    array3D p1 = {{2.0,5.0,7.0}};
    array3D n1 = {{1,0,0}};
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectLinePlane(p0,d0,p1,n1,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectSegmentPlane passed
 */
bool testIntersectSegmentPlane(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{-2,2,-3}};
    array3D d0 = {{4.0/cc,-2.0/cc,3.0/cc}};
    array3D p1 = p0 + 10.0*d0;
    
    array3D p2 = {{2.0,5.0,7.0}};
    array3D n2 = {{1,0,0}};
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectSegmentPlane(p0,p1,p2,n2,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectPlanePlane passed
 */
bool testIntersectPlanePlane(){
    
    bool check(true);
    
    array3D p0 = {{0,0,0}};
    array3D n0 = {{0,1,0}};
    
    array3D p1 = {{0.0,cos(M_PI/3.0),sin(M_PI/3.0)}};
    array3D n1 = {{1,0,0}};
    
    array3D testdir = {{1.0,0.0,0.0}};
    array3D xp, dir;
    
    check &= intersectPlanePlane(p0,n0,p1,n1,xp, dir);
    check &= (dotProduct(xp,n0) < 1.E-12) && (dotProduct(xp,n1) < 1.E-12); //check line origin belongs to planes
    check &= ((std::abs(dotProduct(testdir,dir)) -1.0)< 1.e-12); //test direction, indipendently from sign.
    
    return check;
}

/*!
 * \return true, test intersectPlaneBox passed
 */
bool testIntersectPlaneBox(){

    bool check(true);
    std::vector<array3D> intersects;
    std::vector<array3D> expectedIntersects;

    expectedIntersects.push_back( {{1.0,0.5,0.0}} );
    expectedIntersects.push_back( {{1.0,1.0,0.5}} );
    expectedIntersects.push_back( {{0.5,1.0,0.0}} );

    array3D p0 = {{0.5,1.0,0.0}};
    array3D n0 = {{1.0,1.0,-1.0}};
    n0 /= norm2(n0);

    array3D a0 = {{0.0,0.0,0.0}};
    array3D a1 = {{1.0,1.0,1.0}};

    check = intersectPlaneBox(p0,n0,a0,a1,intersects);
    check &= (intersects.size() == 3);

    if(check){
        for( int i=0; i<3; ++i){
            array3D diff(intersects[i]-expectedIntersects[i]);
            double distance( norm2(diff) );
            check &= utils::DoubleFloatingEqual()( distance, 0 );
        }
    }

    return check;
}

/*!
 * \return true, test intersectLineTriangle passed
 */
bool testIntersectLineTriangle(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{-2,2,-3}};
    array3D d0 = {{4.0/cc,-2.0/cc,3.0/cc}};
    
    array3D p1 = {{2.0,-1.0,-1.0}};
    array3D p2 = {{2.0,2.0,-1.0}};
    array3D p3 = {{2.0,-1.0,1.0}};
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectLineTriangle(p0,d0,p1,p2,p3,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectSegmentTriangle passed
 */
bool testIntersectSegmentTriangle(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{-2,2,-3}};
    array3D d0 = {{4.0/cc,-2.0/cc,3.0/cc}};
    array3D p1 = p0 + 10.0*d0;
    
    array3D p2 = {{2.0,-1.0,-1.0}};
    array3D p3 = {{2.0,2.0,-1.0}};
    array3D p4 = {{2.0,-1.0,1.0}};
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectSegmentTriangle(p0,p1,p2,p3,p4,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectLineSimplex passed
 */
bool testIntersectLineSimplex(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{-2,2,-3}};
    array3D d0 = {{4.0/cc,-2.0/cc,3.0/cc}};
    
    std::vector<array3D> polygon;
    polygon.push_back({{2.0,-1.0,-1.0}});
    polygon.push_back({{2.0,2.0,-1.0}});
    polygon.push_back({{2.0,2.0,1.0}});
    polygon.push_back({{2.0,0.5,2.0}});
    polygon.push_back({{2.0,-1.0,1.0}});
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectLineSimplex(p0,d0,polygon,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectSegmentSimplex passed
 */
bool testIntersectSegmentSimplex(){
    
    bool check(true);
    double cc = pow(29,0.5);
    
    array3D p0 = {{-2,2,-3}};
    array3D d0 = {{4.0/cc,-2.0/cc,3.0/cc}};
    array3D p1 = p0 + 10.0*d0;
    
    std::vector<array3D> polygon;
    polygon.push_back({{2.0,-1.0,-1.0}});
    polygon.push_back({{2.0,2.0,-1.0}});
    polygon.push_back({{2.0,2.0,1.0}});
    polygon.push_back({{2.0,0.5,2.0}});
    polygon.push_back({{2.0,-1.0,1.0}});
    
    array3D test = {{2.0,0.0,0.0}};
    array3D xp;
    
    check &= intersectSegmentSimplex(p0,p1,polygon,xp);
    check &= (norm2(test-xp)<1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectBoxBox passed
 */
bool testIntersectBoxBox(){
    
    bool check(true);
    
    array3D min0 = {{0,0,-1}};
    array3D max0 = {{2,2,5}};
    array3D min1 = {{0.5,0.2,2}};
    array3D max1 = {{4,4,6}};
    
    double res2D = 0.0, res3D=0.0;
    
    array3D min_res, max_res;
    
    check &= intersectBoxBox(min0,max0,min1,max1,min_res, max_res, 2);
    for(int i=0; i<2; ++i){
        res2D += pow((min_res[i] - min1[i]),2);
    }
    check &= (pow(res2D,0.5) < 1.E-12);
    
    check &= intersectBoxBox(min0,max0,min1,max1,min_res, max_res, 3);
    for(int i=0; i<3; ++i){
        res3D += pow((min_res[i] - min1[i]),2);
    }
    check &= (pow(res3D,0.5) < 1.E-12);
    
    return check;
}

/*!
 * \return true, test intersectBoxTriangle passed
 */
bool testIntersectBoxTriangle(){
    
    bool check(true);
    
    array3D min0 = {{-0.5,-0.5,-1}};
    array3D max0 = {{0.5,0.5,1}};
    array3D p1 = {{0,0,0.5}};
    array3D p2 = {{1,0,0.5}};
    array3D p3 = {{0,1,0.5}};

    std::vector<array3D> insects;
    
    check &= intersectBoxTriangle(min0,max0,p1,p2,p3,insects);
    check &= (insects.size() == 1);
    return check;
}

/*!
 * \return true, test intersectSegmentBox passed
 */
bool testIntersectSegmentBox(){
    
    bool check(true);
    
    array3D min0 = {{-0.5,-0.5,-1}};
    array3D max0 = {{0.5,0.5,1}};
    array3D p1 = {{-0.6,0,0.5}};
    array3D p2 = {{0.6,0,0.5}};
    
    std::vector<array3D> insects;

    check &= intersectSegmentBox(p1,p2,min0,max0,3);
    check &= intersectSegmentBox(p1,p2,min0,max0,insects,3);
    check &= (insects.size() == 2);
    
    return check;
}

/*!
 * \return true, test intersectBoxSimplex passed
 */
bool testIntersectBoxSimplex(){
    
    bool check(true);
    
    array3D min0 = {{-0.5,-0.5,-1}};
    array3D max0 = {{0.5,0.5,1}};
    std::vector<array3D> polygon;
    polygon.push_back({{0,0,0.5}});
    polygon.push_back({{0.5,-0.5,0.5}});
    polygon.push_back({{1,0,0.5}});
    polygon.push_back({{1,1,0.5}});
    polygon.push_back({{0,1,0.5}});
    
    std::vector<array3D> insects;
    
    check &= intersectBoxSimplex(min0,max0,polygon,true,true,true,insects,3);
    check &= (insects.size() == 4);
    
    return check;
}

/*!
 * \return true, test intersectPointSegment passed
 */
bool testIntersectPointSegment(){
    
    array3D p0, n0, p1, target1, target2;
    p0={{3,4,5}};
    n0 = {{1.4, -5.6, 9.9}};
    n0 /= norm2(n0);
    p1 = p0 + 5.23*n0;
    target1 = p0 + 1.4453623*n0;
    target2 = p1;
    
    return (intersectPointSegment(target1, p0,p1) && intersectPointSegment(target2, p0,p1) );
}

/*!
 * \return true, test intersectPointTriangle passed
 */
bool testIntersectPointTriangle(){
    
    array3D p0, p1, p2, target1, target2;
    p0 = {{3,4,5}};
    p1 = {{5.5,6,9}};
    p2 = {{-3,0,0}};
    
    target1 = 0.12*p0 + 0.55*p1 + 0.33*p2;
    target2 = 0.5*p0 + 0.0*p1 + 0.5*p2;
    
    return (intersectPointTriangle(target1, p0,p1,p2) && intersectPointTriangle(target2, p0,p1,p2) );
}

/*!
 * \return true, test intersectPointBox passed
 */
bool testIntersectPointBox(){
    
    array3D m0, m1, target1, target2;
    m0 = {{-0.5,-0.5,-0.5}};
    m1 = -1.0*m0;    
    target1 = {{0.23,0.48,-0.35}};
    target2 = {{-0.5,0.2,0.2}};
    
    return (intersectPointBox(target1, m0,m1,3) && intersectPointBox(target2,  m0,m1,3) );
}

/*!
 * \return true, test computeAABBSegment passed
 */
bool testComputeAABBSegment(){
    
    array3D p0,p1,m0, m1, dim, test;
    p0 = {{-2,3,-3}};
    p1 = {{6,10,3}};
    test = {{8,7,6}};
    
    computeAABBSegment(p0,p1,m0,m1);
    dim = m1-m0;
    
    return (norm2(dim-test) < 1.e-12); 
}

/*!
 * \return true, test computeAABBTriangle passed
 */
bool testComputeAABBTriangle(){
    
    array3D p0,p1,p2,m0, m1, dim, test;
    p0 = {{-1,0,0.5}};
    p1 = {{2,0,0.5}};
    p2 = {{0,1.5,0.5}};
    test = {{3,1.5,0}};
    
    computeAABBTriangle(p0,p1,p2,m0,m1);
    dim = m1-m0;
    
    return (norm2(dim-test) < 1.e-12); 
}

/*!
 * \return true, test computeAABBSimplex passed
 */
bool testComputeAABBSimplex(){
    
    std::vector<array3D> polygon;
    array3D m0, m1, dim, test;
    polygon.push_back({{-1,0,0.5}});
    polygon.push_back({{2,0,0.5}});
    polygon.push_back({{2,1.5,0.5}});
    polygon.push_back({{0.23, 1.9567, 0.5}});
    polygon.push_back({{0,1.5,0.5}});
    test = {{3,1.9567,0}};
    
    computeAABBSimplex(polygon,m0,m1);
    dim = m1-m0;
    
    return (norm2(dim-test) < 1.e-12); 
}

/*!
 * \return true, test unionAABB passed
 */
bool testUnionAABB(){
    
    bool check(true);
    
    array3D min0 = {{-1 ,0 ,-2}};
    array3D max0 = {{ 4 ,2 ,5}};
    array3D min1 = {{-2 ,1 ,2}};
    array3D max1 = {{ 3 ,3 ,3}};
    
    array3D test_min,test_max;
    test_min = {{-2,0,-2}};
    test_max = {{4,3,5}};
    
    array3D min_res, max_res;
    
    unionAABB(min0,max0,min1,max1,min_res, max_res);
 
    check &= (norm2(test_min - min_res) < 1.E-12); 
    check &= (norm2(test_max - max_res) < 1.E-12);
    
    std::vector<array3D> mins, maxs;
    
    mins.push_back(min0);
    mins.push_back(min1);
    mins.push_back({{ 5 ,0 ,-1}});
    mins.push_back({{ 8 ,1 ,-4}});

    maxs.push_back(max0);
    maxs.push_back(max1);
    maxs.push_back({{ 7 ,2 ,1}});
    maxs.push_back({{ 10 ,3 ,2}});

    test_min = {{-2,0,-4}};
    test_max = {{10,3,5}};
    
    unionAABB(mins,maxs,min_res, max_res);

    check &= (norm2(test_min - min_res) < 1.E-12);
    check &= (norm2(test_max - max_res) < 1.E-12);
 
    return check;
}

/*!
 * \return true, test intersectionAABB passed
 */
bool testIntersectionAABB(){
    
    bool check(true);
    
    array3D min0 = {{-1 ,0 ,-2}};
    array3D max0 = {{ 4 ,2 ,5}};
    array3D min1 = {{-2 ,1 ,2}};
    array3D max1 = {{ 3 ,3 ,3}};
    
    array3D test_min,test_max;
    test_min = {{-1,1,2}};
    test_max = {{3,2,3}};
    
    array3D min_res, max_res;
    
    intersectionAABB(min0,max0,min1,max1,min_res, max_res);
    
    check &= (norm2(test_min - min_res) < 1.E-12);
    check &= (norm2(test_max - max_res) < 1.E-12);
    return check;
}


/*!
 * \return true, test vertexOfBox passed
 */
bool testVertexOfBox(){

    array3D m0 = {{0,0,0}};
    array3D m1 = {{1,1,1}};
    array3D test2 = {{1,0,0}};

    array3D target;
    int index = 1;
    vertexOfBox(index,m0,m1,target);

    return ( norm2(target-test2) < 1.E-12 );
}

/*!
 * \return true, test edgeOfBox passed
 */
bool testEdgeOfBox(){

    bool check(true);

    array3D m0 = {{0,0,0}};
    array3D m1 = {{1,1,1}};

    array3D test2 = {{0,1,0}};
    array3D test3 = {{1,1,0}};
    
    array3D target2,target3;
    int index = 3;
    
    edgeOfBox(index,m0,m1,target2,target3);
   
    check &= ( norm2(target2-test2) < 1.E-12 );
    check &= ( norm2(target3-test3) < 1.E-12 );

    return check;
}


/*!
 * \return true, test faceOfBox passed
 */
bool testFaceOfBox(){

    bool check(true);

    array3D m0 = {{0,0,0}};
    array3D m1 = {{1,1,1}};
    
    array3D test0 = {{0,0,0}};
    array3D test2 = {{0,1,0}};
    array3D test1 = {{1,0,0}};
    array3D test3 = {{1,1,0}};
    
    array3D target2,target3, target0, target1;
    int index = 4;
    
    faceOfBox(index,m0,m1,target0,target1, target3, target2);
    
    check &= ( norm2(target0-test0) < 1.E-12 );
    check &= ( norm2(target1-test1) < 1.E-12 );
    check &= ( norm2(target2-test2) < 1.E-12 );
    check &= ( norm2(target3-test3) < 1.E-12 );
    
    return check;
    
}

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(int argc, char *argv[])
{
    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variabels
    bool                            pass = false;

    // ====================================================================== //
    // RUN SUB-TESTS                                                          //
    // ====================================================================== //
    try {
        std::cout << "Testing convertBarycentricToFlagSegment...";
        pass = testConvertBarycentricToFlagSegment();
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testConvertBarycentricToFlagTriangle();
        std::cout << "Testing convertBarycentricToFlagTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testConvertBarycentricToFlagSimplex();
        std::cout << "Testing convertBarycentricToFlagSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testReconstructPointFromBarycentricSegment();
        std::cout << "Testing reconstructPointFromBarycentricSegment...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testReconstructPointFromBarycentricTriangle();
        std::cout << "Testing reconstructPointFromBarycentricTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testReconstructPointFromBarycentricSimplex();
        std::cout << "Testing reconstructPointFromBarycentricSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testProjectPointLine();
        std::cout << "Testing projectPointLine...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testProjectPointPlane();
        std::cout << "Testing projectPointPlane...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testProjectPointSegment();
        std::cout << "Testing projectPointSegment...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testProjectPointTriangle();
        std::cout << "Testing projectPointTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testRestrictPointTriangle();
        std::cout << "Testing restrictPointTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testProjectCloudTriangle();
        std::cout << "Testing projectCloudTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testProjectPointSimplex();
        std::cout << "Testing projectPointSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testProjectPointCone();
        std::cout << "Testing projectPointCone...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistancePointLine();
        std::cout << "Testing distancePointLine...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistancePointPlane();
        std::cout << "Testing distancePointPlane...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistancePointSegment();
        std::cout << "Testing distancePointSegment...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistancePointTriangle();
        std::cout << "Testing distancePointTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistancePointCone();
        std::cout << "Testing distancePointCone...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistanceCloudTriangle();
        std::cout << "Testing distanceCloudTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistancePointSimplex();
        std::cout << "Testing distancePointSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistanceCloudSimplex();
        std::cout << "Testing distanceCloudSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testDistanceLineLine();
        std::cout << "Testing distanceLineLine...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectLineLine();
        std::cout << "Testing intersectLineLine...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectSegmentSegment();
        std::cout << "Testing intersectSegmentSegment...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectLinePlane();
        std::cout << "Testing intersectLinePlane...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectSegmentPlane();
        std::cout << "Testing intersectSegmentPlane...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectPlanePlane();
        std::cout << "Testing intersectPlanePlane...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectPlaneBox();
        std::cout << "Testing intersectPlaneBox...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectLineTriangle();
        std::cout << "Testing intersectLineTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectSegmentTriangle();
        std::cout << "Testing intersectSegmentTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectLineSimplex();
        std::cout << "Testing intersectLineSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectSegmentSimplex();
        std::cout << "Testing intersectSegmentSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectBoxBox();
        std::cout << "Testing intersectBoxBox...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectBoxTriangle();
        std::cout << "Testing intersectBoxTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectSegmentBox();
        std::cout << "Testing intersectSegmentBox...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectBoxSimplex();
        std::cout << "Testing intersectBoxSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectPointSegment();
        std::cout << "Testing intersectPointSegment...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectPointTriangle();
        std::cout << "Testing intersectPointTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectPointBox();
        std::cout << "Testing intersectPointBox...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testComputeAABBSegment();
        std::cout << "Testing computeAABBSegment...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testComputeAABBTriangle();
        std::cout << "Testing computeAABBTriangle...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testComputeAABBSimplex();
        std::cout << "Testing computeAABBSimplex...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testUnionAABB();
        std::cout << "Testing unionAABB...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testIntersectionAABB();
        std::cout << "Testing intersectionAABB...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testVertexOfBox();
        std::cout << "Testing vertexOfBox...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testEdgeOfBox();
        std::cout << "Testing edgeOfBox...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

        pass = testFaceOfBox();
        std::cout << "Testing faceOfBox...";
        if (!pass) {
            std::cout << " Failed" << std::endl;
            return 1;
        } else {
            std::cout << " Passed" << std::endl;
        }

    } catch (const std::exception &exception) {
        std::cout << exception.what();
        exit(1);
    }

    // ====================================================================== //
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return 0;
}
