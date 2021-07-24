#include <sstream>
#include <sys/stat.h>
#include "Node.hpp"
#include "Element.hpp"
#include "Recorder.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Paraview Overload constructor.
Recorder::Recorder(std::string file, std::string name, unsigned int nparaview, unsigned int nsample, unsigned int precision) :
File(file), Name(name), nSample(nsample), Counter(0), nParaview(nparaview), Precision(precision) { 
    //Does nothing.
}

//Point/Element Overload constructor.
Recorder::Recorder(std::string file, std::string name, std::string type, std::vector<unsigned int> index, unsigned int nsample, unsigned int precision) :
File(file), Name(name), Response(type), nSample(nsample), Counter(0), Precision(precision), IDs(index){ 
    //Does nothing.
}

//Section Overload constructor.
Recorder::Recorder(std::string file, std::string name, std::string type, std::vector<double> coordinates, std::vector<unsigned int> index, unsigned int nsample, unsigned int precision) : 
File(file), Name(name), Response(type), nSample(nsample), Counter(0), Precision(precision), Position(coordinates), IDs(index){ 
    //In case the section is an area.
    if(Position.size() == 1)
        Position.push_back(0.0);
}

//Virtual destructor.
Recorder::~Recorder(){
    //Does nothing.
}

//Clone the 'recorder' object.
std::unique_ptr<Recorder>
Recorder::CopyRecorder(){
    //Selects the recorder type.
    if(strcasecmp(Name.c_str(),"PARAVIEW") == 0){
        return std::make_unique<Recorder>(File, Name, nParaview, nSample, Precision);
    }
    else if(strcasecmp(Name.c_str(),"SECTION") == 0){
        return std::make_unique<Recorder>(File, Name, Response, Position, IDs, nSample, Precision);
    }
    else{
        return std::make_unique<Recorder>(File, Name, Response, IDs, nSample, Precision);
    }
}

//Return the recorder name
std::string 
Recorder::GetName(){
    return Name;
}

//Initialize the recorder.
void 
Recorder::Initialize(std::shared_ptr<Mesh> &mesh, unsigned int nsteps){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    if(strcasecmp(Name.c_str(),"PARAVIEW") == 0){
        //Gets node and element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

        //Local node indexes mapping.
        unsigned int k = 0;
        for(auto it : Nodes){
            auto &nTag = it.first;
            Tag[nTag]  = k;
            k++;
        }
    }
    else if(strcasecmp(Name.c_str(),"NODE") == 0){       
        //Update the file path if is element recorder
        std::string secondPath = GetSpacedName(filePath, " ");
        std::string theFile = secondPath + "/../Solution/" + Combo + "/" + File;

        //Creates the output file.
        OutputFile.open(theFile.c_str()); 
        OutputFile.precision(Precision);
        OutputFile.setf(std::ios::scientific);

        //Gets node information from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

        //Number of degree-of-freedom for this partition.
        unsigned int nTotalDofsPartition = 0;
        for(unsigned int k = 0; k < IDs.size(); k++)
            nTotalDofsPartition += Nodes[IDs[k]]->GetNumberOfDegreeOfFreedom();

        //Number of nodes in this partition.
        OutputFile << IDs.size() << " " << nTotalDofsPartition << " " << numberOfTotalDofs << " " << nsteps << "\n";

        //List with node information.
        for(unsigned int k = 0; k < IDs.size(); k++){
            OutputFile << IDs[k] << " " << Nodes[IDs[k]]->GetNumberOfDegreeOfFreedom();

            std::vector<int> TotalDofs = Nodes[IDs[k]]->GetTotalDegreeOfFreedom();
            for(unsigned int i = 0; i < TotalDofs.size(); i++){
                OutputFile << " " << TotalDofs[i];
            }

            OutputFile << "\n";
        }

    }
    else if(strcasecmp(Name.c_str(),"ELEMENT") == 0){
        //Update the file path if is element recorder
        std::string secondPath = GetSpacedName(filePath, " ");
        std::string theFile = secondPath + "/../Solution/" + Combo + "/" + File;

        //Creates the output file.
        OutputFile.open(theFile.c_str()); 
        OutputFile.precision(Precision);
        OutputFile.setf(std::ios::scientific);

        //Gets node and element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

        //Gets node and element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

        //Number of nodes in this partition.
        OutputFile << IDs.size() << " " << numberOfTotalDofs << " " << nsteps << "\n";

        //List with element information.
        if(strcasecmp(Response.c_str(),"INTERNALFORCE") == 0){
                
            for(unsigned int k = 0; k < IDs.size(); k++){
                OutputFile << IDs[k] << " " << Elements[IDs[k]]->GetNumberOfDegreeOfFreedom();

                //Returns the Node Connectivity Indeces.
                std::vector<unsigned int> conn = Elements[IDs[k]]->GetNodes();

                for(unsigned int j = 0; j < conn.size(); j++){
                    std::vector<int> TotalDofs = Nodes[conn[j]]->GetTotalDegreeOfFreedom();

                    for(unsigned int i = 0; i < TotalDofs.size(); i++)
                        OutputFile << " " << TotalDofs[i]; 
                }

                OutputFile << "\n";
            }
        }
        else if(strcasecmp(Response.c_str(),"STRESS") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd MatInfo = Elements[IDs[k]]->GetStress();

                OutputFile << IDs[k] << " " << MatInfo.rows() << " " << MatInfo.cols() << "\n";
            }
        }
        else if(strcasecmp(Response.c_str(),"STRAIN") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd MatInfo = Elements[IDs[k]]->GetStrain();

                OutputFile << IDs[k] << " " << MatInfo.rows() << " " << MatInfo.cols() << "\n";
            }
        }
        else if(strcasecmp(Response.c_str(),"STRAINRATE") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd MatInfo = Elements[IDs[k]]->GetStrainRate();

                OutputFile << IDs[k] << " " << MatInfo.rows() << " " << MatInfo.cols() << "\n";
            }
        }
    }
    else if(strcasecmp(Name.c_str(),"SECTION") == 0){
        //Update the file path if is element recorder
        std::string secondPath = GetSpacedName(filePath, " ");
        std::string theFile = secondPath + "/../Solution/" + Combo + "/" + File;

        //Creates the output file.
        OutputFile.open(theFile.c_str()); 
        OutputFile.precision(Precision);
        OutputFile.setf(std::ios::scientific);

        //Gets node and element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

        //Number of section in this partition.
        OutputFile << IDs.size() << " " << nsteps;
        for(unsigned int k = 0; k < Position.size(); k++)
            OutputFile << " " << Position[k];
        OutputFile << "\n";

        //TODO: Section recorder.
        if(strcasecmp(Response.c_str(),"STRAIN") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd SectionInfo = Elements[IDs[k]]->GetStrainAt();
                OutputFile << IDs[k] << " " << SectionInfo.rows() << " " << SectionInfo.cols() << "\n";
            }
        }
        else if(strcasecmp(Response.c_str(),"STRESS") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd SectionInfo = Elements[IDs[k]]->GetStressAt();
                OutputFile << IDs[k] << " " << SectionInfo.rows() << " " << SectionInfo.cols() << "\n";
            }
        }
    }
}

//Writes state variable data to the file.
void 
Recorder::WriteResponse(std::shared_ptr<Mesh> &mesh, unsigned int step){
    //Starts profiling this fuction.
    PROFILE_FUNCTION();

    //Selects the object to be written.
    if(strcasecmp(Name.c_str(),"PARAVIEW") == 0){
        WriteVTKFiles(mesh, step);
    }
    else if(strcasecmp(Name.c_str(),"NODE") == 0){
        WriteNodalResponse(mesh);
    }
    else if(strcasecmp(Name.c_str(),"ELEMENT") == 0){
        WriteElementResponse(mesh);
    }
    else if(strcasecmp(Name.c_str(),"SECTION") == 0){
        WriteSectionResponse(mesh);
    }
}

//Close the recorder.
void 
Recorder::Finalize(){
    //Closes the output File:
    if(strcasecmp(Name.c_str(),"PARAVIEW") != 0)
        OutputFile.close(); 
}

//Sets the combination's name.
void 
Recorder::SetComboName(std::string name){
    Combo = name;
}

//Writes state variable data to the file.
void 
Recorder::WriteNodalResponse(std::shared_ptr<Mesh> &mesh){
    //Gets all nodes from the mesh.
    std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

    //Some nodes are to be recorder: all degree-of-freedom are written.
    for(unsigned int k = 0; k < IDs.size(); k++){
        //Vector that stores solution.
        Eigen::VectorXd U;

        //Selects the state variable according to this recorder.
        if(strcasecmp(Response.c_str(),"DISP") == 0){
            U = Nodes[IDs[k]]->GetDisplacements();
        }
        else if(strcasecmp(Response.c_str(),"VEL") == 0){
            U = Nodes[IDs[k]]->GetVelocities();
        }
        else if(strcasecmp(Response.c_str(),"ACCEL") == 0){
            U = Nodes[IDs[k]]->GetAccelerations();
        }
        else if(strcasecmp(Response.c_str(),"REACTION") == 0){
            U = Nodes[IDs[k]]->GetReaction();
        }

        //Writes the node state.
        for(unsigned int j = 0; j < U.size(); j++)
            OutputFile << U(j) << " ";
    }

    OutputFile << "\n"; 
}

//Writes state variable data to the file.
void 
Recorder::WriteElementResponse(std::shared_ptr<Mesh> &mesh){
    //Gets node and element information from the mesh.
    std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

    //Select the element response to be recorded.
    if(strcasecmp(Response.c_str(),"STRAIN") == 0){
        for(unsigned int k = 0; k < IDs.size(); k++){
            //Gets the material strain.
            Eigen::MatrixXd Strain = Elements[IDs[k]]->GetStrain();

            //Writes the internal force vector.
            for(unsigned int i = 0; i < Strain.rows(); i++){
                for(unsigned int j = 0; j < Strain.cols(); j++)
                    OutputFile << Strain(i,j) << " ";
            }
        }
    }
    else if(strcasecmp(Response.c_str(),"STRESS") == 0){
        for(unsigned int k = 0; k < IDs.size(); k++){
            //Gets the material stress.
            Eigen::MatrixXd Stress = Elements[IDs[k]]->GetStress();

            //Writes the internal force vector.
            for(unsigned int i = 0; i < Stress.rows(); i++){
                for(unsigned int j = 0; j < Stress.cols(); j++)
                    OutputFile << Stress(i,j) << " ";
            }
        }
    }
    else if(strcasecmp(Response.c_str(),"INTERNALFORCE") == 0){
        for(unsigned int k = 0; k < IDs.size(); k++){
            //Internal force vector.
            Eigen::VectorXd force = Elements[IDs[k]]->ComputeInternalForces();

            //Writes the internal force vector.
            for(unsigned int i = 0; i < force.size(); i++)
                OutputFile << force(i) << " ";
        }
    }
    else if(strcasecmp(Response.c_str(),"STRAINRATE") == 0){
        for(unsigned int k = 0; k < IDs.size(); k++){
            //Gets the material stress.
            Eigen::MatrixXd StrainRate = Elements[IDs[k]]->GetStrainRate();

            //Writes the internal force vector.
            for(unsigned int i = 0; i < StrainRate.rows(); i++){
                for(unsigned int j = 0; j < StrainRate.cols(); j++)
                    OutputFile << StrainRate(i,j) << " ";
            }
        }
    }
    
    OutputFile << "\n"; 
}

//Writes section element data to the file.
void 
Recorder::WriteSectionResponse(std::shared_ptr<Mesh> &mesh){
    //Gets all elements from the mesh.
    std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

    //Select the element response to be recorded.
    if(strcasecmp(Response.c_str(),"STRAIN") == 0){
        for(unsigned int k = 0; k < IDs.size(); k++){
            //Gets the section strain matrix.
            Eigen::MatrixXd SecStrain = Elements[IDs[k]]->GetStrainAt(Position[0], Position[1]);

            //Writes the strain vector at each integration section and position.
            for(unsigned int i = 0; i < SecStrain.rows(); i++){
                for(unsigned int j = 0; j < SecStrain.cols(); j++)
                    OutputFile << SecStrain(i,j) << " ";
            }
        }
    }
    else if(strcasecmp(Response.c_str(),"STRESS") == 0){
        for(unsigned int k = 0; k < IDs.size(); k++){
            //Gets the section strain matrix.
            Eigen::MatrixXd SecStress = Elements[IDs[k]]->GetStressAt(Position[0], Position[1]);

            //Writes the strain vector at each integration section and position.
            for(unsigned int i = 0; i < SecStress.rows(); i++){
                for(unsigned int j = 0; j < SecStress.cols(); j++)
                    OutputFile << SecStress(i,j) << " ";
            }
        }
    }
    OutputFile << "\n"; 
}

//Writes output data in VTK format to the file
void 
Recorder::WriteVTKFiles(std::shared_ptr<Mesh> &mesh, unsigned int step){
    if(step % nSample == 0){
        //Creates the paraview output file
        std::stringstream iter;
        iter << Counter;

        //Header line skips for XML file format
        std::string head1("  ");
        std::string head2("    ");
        std::string head3("      ");
        std::string head4("        ");

        //Update the file path for the paraview output
        std::string firstPath = GetSpacedName(filePath, " ");
        std::string Paraview  = firstPath  + "/../Paraview/" +  Combo + "/" + File + "." + iter.str() + ".vtu";

        //Creates the output file
        OutputFile.open(Paraview.c_str()); 
        OutputFile.precision(Precision);
        OutputFile.setf(std::ios::scientific);

        //Gets all Points/Elements from the mesh
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();
        std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

        //VTK header file information
        OutputFile << "<?xml version=\"1.0\"?>\n";
        OutputFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"BigEndian\">\n" << head1;
        OutputFile << "<UnstructuredGrid>\n" << head2 ;
        OutputFile << "<Piece NumberOfPoints=\"" << Nodes.size() << "\" NumberOfCells=\"" << Elements.size() << "\">\n" << head3;

        //Node Coordinates
        OutputFile << "<Points>\n" << head4;
        OutputFile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
        
        for(auto it : Nodes){
            auto &nTag = it.first;

            //Get Node information
            Eigen::VectorXd X = Nodes[nTag]->GetCoordinates();

            if(nDimensions == 2)
                OutputFile << head4 << X(0) << " " << X(1) << " 0.00000\n";
            else if(nDimensions == 3)
                OutputFile << head4 << X(0) << " " << X(1) << " " << X(2) << "\n";
        }

        OutputFile << head4; 
        OutputFile << "</DataArray>\n" << head3;
        OutputFile << "</Points>\n" << head3;

        //Elements Connectivities
        OutputFile << "<Cells>\n" << head4;
        OutputFile << "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";

        for(auto it : Elements){
            auto &eTag = it.first;
            std::vector<unsigned int> connection = Elements[eTag]->GetNodes();

            OutputFile << head4;
            for(unsigned int k = 0; k < connection.size(); k++)
                OutputFile << Tag[connection[k]] << " ";
            OutputFile << "\n";
        }

        OutputFile << head4;
        OutputFile << "</DataArray>\n" << head4;

        //Elements Number of Nodes Offset.
        OutputFile << "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n" << head4;

        unsigned int nOffsets = 0;
        for(auto it : Elements){
            auto &eTag = it.first;

            nOffsets += Elements[eTag]->GetNumberOfNodes();
            OutputFile << " " << nOffsets;
        }

        OutputFile << "\n" << head4;
        OutputFile << "</DataArray>\n" << head4;

        //Elements Geometry Types.
        OutputFile << "<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n" << head4;

        for(auto it : Elements){
            auto &eTag = it.first;
            OutputFile << Elements[eTag]->GetVTKCellType() << " ";
        }

        OutputFile << "\n" << head4;
        OutputFile << "</DataArray>\n" << head3;
        OutputFile << "</Cells>\n" << head3;

        //Writes the Node information
        OutputFile << "<PointData>\n" << head4;

        //Writes the node identifiers
        OutputFile << "<DataArray type=\"Int32\" Name=\"GlobalNodeId\" format=\"ascii\">\n" << head4;

        for(auto it : Nodes){
            auto &nTag = it.first;
            OutputFile << nTag << " ";
        }
        
        OutputFile << "\n" << head4;
        OutputFile << "</DataArray>\n" << head4;
        
        //Writes the displacement at each node.
        OutputFile << "<DataArray type=\"Float32\" Name=\"Displacements\" NumberOfComponents=\"3\" ComponentName0=\"Ux\" ComponentName1=\"Uy\" ComponentName2=\"Uz\" Format=\"ascii\">\n";

        for(auto it : Nodes){
            auto &nTag = it.first;
            //Get Node displacements.
            Eigen::VectorXd U = Nodes[nTag]->GetDisplacements();
            SetThreshold(U);

            if(nDimensions == 2)
                OutputFile << head4 << U(0) << " " << U(1) << " 0.00000\n";
            else if(nDimensions == 3)
                OutputFile << head4 << U(0) << " " << U(1) << " " << U(2) << "\n";
        }

        OutputFile << head4;
        OutputFile << "</DataArray>\n" << head4;

        //Writes the velocity at each node.
        OutputFile << "<DataArray type=\"Float32\" Name=\"Velocities\" NumberOfComponents=\"3\" ComponentName0=\"Vx\" ComponentName1=\"Vy\" ComponentName2=\"Vz\" Format=\"ascii\">\n";

        for(auto it : Nodes){
            auto &nTag = it.first;
            //Get Node velocities.
            Eigen::VectorXd V = Nodes[nTag]->GetVelocities();
            SetThreshold(V);

            if(nDimensions == 2)
                OutputFile << head4 << V(0) << " " << V(1) << " 0.00000\n";
            else if(nDimensions == 3)
                OutputFile << head4 << V(0) << " " << V(1) << " " << V(2) << "\n";
        }

        OutputFile << head4;
        OutputFile << "</DataArray>\n" << head4;

        //Writes the acceleration at each node.
        OutputFile << "<DataArray type=\"Float32\" Name=\"Accelerations\" NumberOfComponents=\"3\" ComponentName0=\"Ax\" ComponentName1=\"Ay\" ComponentName2=\"Az\" Format=\"ascii\">\n";
        
        for(auto it : Nodes){
            auto &nTag = it.first;
            //Get Node accelerations
            Eigen::VectorXd A = Nodes[nTag]->GetAccelerations();
            SetThreshold(A);

            if(nDimensions == 2)
                OutputFile << head4 << A(0) << " " << A(1) << " 0.00000\n";
            else if(nDimensions == 3)
                OutputFile << head4 << A(0) << " " << A(1) << " " << A(2) << "\n";
        }

        OutputFile << head4;
        OutputFile << "</DataArray>\n" << head3;
        OutputFile << "</PointData>\n" << head3;

        //Writes the Element information
        OutputFile << "<CellData>\n" << head4;

        //Writes the element identifiers
        OutputFile << "<DataArray type=\"Int32\" Name=\"GlobalElementId\" format=\"ascii\">\n";

        for(auto it : Elements){
            auto &eTag = it.first;
            OutputFile << eTag << " ";
        }
        
        OutputFile << "\n" << head4;
        OutputFile << "</DataArray>\n" << head4;

        //Writes the element group identifier
        OutputFile << "<DataArray type=\"Int32\" Name=\"GroupElementId\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"100\">\n";

        for(auto it : Elements){
            auto &eTag = it.first;
            OutputFile << Elements[eTag]->GetVTKGroupType() << " ";
        }
        
        OutputFile << "\n" << head4;
        OutputFile << "</DataArray>\n" << head4;

        //Writes the strains at each element
        OutputFile << "<DataArray type=\"Float32\" Name=\"ElementStrains\" NumberOfComponents=\"18\" ComponentName0=\"E11 (Solid)\" ComponentName1=\"E22 (Solid)\" ComponentName2=\"E33 (Solid)\" ComponentName3=\"E12 (Solid)\" ComponentName4=\"E23 (Solid)\" ComponentName5=\"E13 (Solid)\" ComponentName6=\"e11 (Shell)\" ComponentName7=\"e22 (Shell)\" ComponentName8=\"e12 (Shell)\" ComponentName9=\"kappa11 (Shell)\" ComponentName10=\"kappa22 (Shell)\" ComponentName11=\"kappa12 (Shell)\" ComponentName12=\"epsilon (Frame)\" ComponentName13=\"gamma2 (Frame)\" ComponentName14=\"gamma3 (Frame)\" ComponentName15=\"phi (Frame)\" ComponentName16=\"kappa2 (Frame)\" ComponentName17=\"kappa3 (Frame)\" Format=\"ascii\">\n";

        for(auto it : Elements){
            auto &eTag = it.first;

            //Gets the material stress.
            Eigen::VectorXd Strain = Elements[eTag]->GetVTKResponse("Strain");
            SetThreshold(Strain);

            //PATCH: Sets the DRM strains to be zero (visualization)
            if(IsDRMElem[eTag]){
                Strain.fill(0.0);
            }

            OutputFile << head4 << Strain.transpose() << "\n";
        }

        OutputFile << head4;
        OutputFile << "</DataArray>\n" << head4;

        //Writes the stresses at each element
        OutputFile << "<DataArray type=\"Float32\" Name=\"ElementStresses\" NumberOfComponents=\"18\" ComponentName0=\"S11 (Solid)\" ComponentName1=\"S22 (Solid)\" ComponentName2=\"S33 (Solid)\" ComponentName3=\"S12 (Solid)\" ComponentName4=\"S23 (Solid)\" ComponentName5=\"S13 (Solid)\" ComponentName6=\"N11 (Shell)\" ComponentName7=\"N22 (Shell)\" ComponentName8=\"N12 (Shell)\" ComponentName9=\"M11 (Shell)\" ComponentName10=\"M22 (Shell)\" ComponentName11=\"M12 (Shell)\" ComponentName12=\"N (Frame)\" ComponentName13=\"V2 (Frame)\" ComponentName14=\"V3 (Frame)\" ComponentName15=\"T (Frame)\" ComponentName16=\"M2 (Frame)\" ComponentName17=\"M3 (Frame)\" Format=\"ascii\">\n";

        for(auto it : Elements){
            auto &eTag = it.first;

            //Gets the material stress.
            Eigen::VectorXd Stress = Elements[eTag]->GetVTKResponse("Stress");
            SetThreshold(Stress);

            //PATCH: Sets the DRM strains to be zero (visualization)
            if(IsDRMElem[eTag]){
                Stress.fill(0.0);
            }

            OutputFile << head4 << Stress.transpose() << "\n";
        }

        OutputFile << head4;
        OutputFile << "</DataArray>\n" << head3;
        
        OutputFile << "</CellData>\n"  << head2;

        //VTK footer file information:
        OutputFile << "</Piece>\n" << head1;
        OutputFile << "</UnstructuredGrid>\n";
        OutputFile << "</VTKFile>";

        OutputFile.close();
        Counter++;
    }
}

//Fix blank spaces provided by user in path.
std::string 
Recorder::GetSpacedName(std::string theFile, std::string toReplace){
    //Auxiliary variables.
    std::string auxName = theFile;
    std::string correctedFile;

    //Replace strings.
    size_t pos;
    pos = auxName.find("~");
    while(pos != std::string::npos){
        auxName.replace(pos, std::string("~").length(), toReplace);
        pos = auxName.find("~");
    }

    //Modified File Name.
    correctedFile = auxName;

    return correctedFile;
}

//Sets the possible DRM Element indexes. 
void 
Recorder::SetDRMParaviewInterface(std::map<unsigned int, bool>& DRMElems){
    IsDRMElem = DRMElems;
}

//Windows patch for small values
void 
Recorder::SetThreshold(Eigen::VectorXd &U, double WinTol){
    for(unsigned int k = 0; k < U.size(); k++){
        if(fabs(U(k)) < WinTol)
            U(k) = 0.0;
    }
}
