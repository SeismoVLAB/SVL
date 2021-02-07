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
    if(strcasecmp(Name.c_str(),"paraview") == 0){
        return std::make_unique<Recorder>(File, Name, nParaview, nSample, Precision);
    }
    else if(strcasecmp(Name.c_str(),"section") == 0){
        return std::make_unique<Recorder>(File, Name, Response, Position, IDs, nSample, Precision);
    }
    else{
        return std::make_unique<Recorder>(File, Name, Response, IDs, nSample, Precision);
    }
}

//Initialize the recorder.
void 
Recorder::Initialize(std::shared_ptr<Mesh> &mesh, unsigned int nsteps){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    if(strcasecmp(Name.c_str(),"paraview") == 0){
        //Creates the folder where solutions are stored.
        MakeFolder("/../Paraview/");

        //Gets node and element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

        //Local node indeces mapping.
        unsigned int k = 0;
        for(auto it : Nodes){
            auto &nTag = it.first;
            Tag[nTag]  = k;
            k++;
        }
    }
    else if(strcasecmp(Name.c_str(),"node") == 0){
        //Creates the folder where solutions are stored.
        MakeFolder("/../Solution/");
        
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
    else if(strcasecmp(Name.c_str(),"element") == 0){
        //Creates the folder where solutions are stored.
        MakeFolder("/../Solution/");

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
        if(strcasecmp(Response.c_str(),"internalforce") == 0){
                
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
        else if(strcasecmp(Response.c_str(),"stress") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd MatInfo = Elements[IDs[k]]->GetStress();

                OutputFile << IDs[k] << " " << MatInfo.rows() << " " << MatInfo.cols() << "\n";
            }
        }
        else if(strcasecmp(Response.c_str(),"strain") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd MatInfo = Elements[IDs[k]]->GetStrain();

                OutputFile << IDs[k] << " " << MatInfo.rows() << " " << MatInfo.cols() << "\n";
            }
        }
        else if(strcasecmp(Response.c_str(),"strainrate") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd MatInfo = Elements[IDs[k]]->GetStrainRate();

                OutputFile << IDs[k] << " " << MatInfo.rows() << " " << MatInfo.cols() << "\n";
            }
        }
    }
    else if(strcasecmp(Name.c_str(),"section") == 0){
        //Creates the folder where solutions are stored.
        MakeFolder("/../Solution/");

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
        OutputFile << IDs.size() << " " << nsteps << " " << Position[0] << " " << Position[1] << "\n";

        //TODO: Section recorder.
        if(strcasecmp(Response.c_str(),"strain") == 0){
            for(unsigned int k = 0; k < IDs.size(); k++){
                Eigen::MatrixXd SectionInfo = Elements[IDs[k]]->GetStrainAt();
                OutputFile << IDs[k] << " " << SectionInfo.rows() << " " << SectionInfo.cols() << "\n";
            }
        }
        else if(strcasecmp(Response.c_str(),"stress") == 0){
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
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Selects the object to be written.
    if(strcasecmp(Name.c_str(),"paraview") == 0){
        WriteVTKFiles(mesh, step);
    }
    else if(strcasecmp(Name.c_str(),"node") == 0){
        WriteNodalResponse(mesh);
    }
    else if(strcasecmp(Name.c_str(),"element") == 0){
        WriteElementResponse(mesh);
    }
    else if(strcasecmp(Name.c_str(),"section") == 0){
        WriteSectionResponse(mesh);
    }
}

//Close the recorder.
void 
Recorder::Finalize(){
    //Closes the output File:
    if(strcasecmp(Name.c_str(),"paraview") != 0)
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
        if(strcasecmp(Response.c_str(),"disp") == 0){
            U = Nodes[IDs[k]]->GetDisplacements();
        }
        else if(strcasecmp(Response.c_str(),"vel") == 0){
            U = Nodes[IDs[k]]->GetVelocities();
        }
        else if(strcasecmp(Response.c_str(),"accel") == 0){
            U = Nodes[IDs[k]]->GetAccelerations();
        }
        else if(strcasecmp(Response.c_str(),"reaction") == 0){
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
    if(strcasecmp(Response.c_str(),"strain") == 0){
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
    else if(strcasecmp(Response.c_str(),"stress") == 0){
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
    else if(strcasecmp(Response.c_str(),"internalforce") == 0){
        for(unsigned int k = 0; k < IDs.size(); k++){
            //Internal force vector.
            Eigen::VectorXd force = Elements[IDs[k]]->ComputeInternalForces();

            //Writes the internal force vector.
            for(unsigned int i = 0; i < force.size(); i++)
                OutputFile << force(i) << " ";
        }
    }
    else if(strcasecmp(Response.c_str(),"strainrate") == 0){
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
    if(strcasecmp(Response.c_str(),"strain") == 0){
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
    else if(strcasecmp(Response.c_str(),"stress") == 0){
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

//Writes output data in VTK format to the file.
void 
Recorder::WriteVTKFiles(std::shared_ptr<Mesh> &mesh, unsigned int step){
    if((step-1) % nSample == 0){
        //Creates the paraview output file.
        std::stringstream iter;
        iter << Counter;

        //Update the file path for the paraview output
        std::string firstPath = GetSpacedName(filePath, " ");
        std::string Paraview  = firstPath  + "/../Paraview/" +  Combo + "/" + File + "." + iter.str() + ".vtk";

        //Creates the output file.
        OutputFile.open(Paraview.c_str()); 
        OutputFile.precision(Precision);
        OutputFile.setf(std::ios::scientific);

        //Model Global Information:
        OutputFile << "# vtk DataFile Version 4.0\n";
        OutputFile << "SEISMOVLAB: VTK POST-PROCESS FILE\n";
        OutputFile << "ASCII\n";
        OutputFile << "DATASET UNSTRUCTURED_GRID\n";
        OutputFile << "\n";

        //Gets all Points/Elements from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();
        std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

        //Animated nodal mesh coordinates.
        OutputFile << "POINTS " << Nodes.size() << " float\n";
        for(auto it : Nodes){
            auto &nTag = it.first;

            //Get Node information.
            Eigen::VectorXd X = Nodes[nTag]->GetCoordinates();

            OutputFile << X(0) << " " << X(1) << " ";
            if(X.size() == 2)
                OutputFile << "0.00000\n";
            else if(X.size() == 3)
                OutputFile << X(2) << "\n";
        }

        //Elements Connectivities.
        OutputFile << "\nCELLS " << Elements.size() << " " << nParaview << "\n";
        for(auto it : Elements){
            auto &eTag = it.first;
            std::vector<unsigned int> connection = Elements[eTag]->GetNodes();

            OutputFile << connection.size();
            for(unsigned int k = 0; k < connection.size(); k++)
                OutputFile << " " << Tag[connection[k]];
            OutputFile << "\n";
        }

        //Elements Geometry Types.
        OutputFile << "\nCELL_TYPES " << Elements.size() << "\n";
        for(auto it : Elements){
            auto &eTag = it.first;
            OutputFile << Elements[eTag]->GetVTKCellType() <<"\n";
        }

        //Writes the displacement at each node.
        OutputFile << "\nPOINT_DATA " << Nodes.size() << "\n";
        OutputFile << "VECTORS Displacement float\n";

        for(auto it : Nodes){
            auto &nTag = it.first;
            //Get Node information.
            Eigen::VectorXd U = Nodes[nTag]->GetDisplacements();

            OutputFile << U(0) << " " << U(1) << " ";
            if(nDimensions == 2)
                OutputFile << "0.00000\n";
            else if(nDimensions == 3)
                OutputFile << U(2) << "\n";
        }

        //Writes the velocity at each node.
        OutputFile << "\nVECTORS Velocity float\n";

        for(auto it : Nodes){
            auto &nTag = it.first;

            //Get Node information.
            Eigen::VectorXd V = Nodes[nTag]->GetVelocities();

            OutputFile << V(0) << " " << V(1) << " ";
            if(nDimensions == 2)
                OutputFile << "0.00000\n";
            else if(nDimensions == 3)
                OutputFile << V(2) << "\n";
        }

        //Writes the acceleration at each node.
        OutputFile << "\nVECTORS Acceleration float\n";

        for(auto it : Nodes){
            auto &nTag = it.first;

            //Get Node information.
            Eigen::VectorXd A = Nodes[nTag]->GetAccelerations();

            OutputFile << A(0) << " " << A(1) << " ";
            if(nDimensions == 2)
                OutputFile << "0.00000\n";
            else if(nDimensions == 3)
                OutputFile << A(2) << "\n";
        }

        //Writes the acceleration at each node.
        OutputFile << "\nVECTORS Reaction float\n";

        for(auto it : Nodes){
            auto &nTag = it.first;
            //Get Node information.
            Eigen::VectorXd R = Nodes[nTag]->GetReaction();

            OutputFile << R(0) << " " << R(1) << " ";
            if(nDimensions == 2)
                OutputFile << "0.00000\n";
            else if(nDimensions == 3)
                OutputFile << R(2) << "\n";
        }

        //TODO: Implement Paraview Record for Stresses/Internal Forces.
        OutputFile << "\nCELL_DATA " << Elements.size() << "\n"; 
        OutputFile << "FIELD attributes 2\n";

        //Writes the acceleration at each node.
        OutputFile << "Strains 6 " << Elements.size() << " float\n";
        for(auto it : Elements){
            auto &eTag = it.first;

            //Gets the material stress.
            Eigen::VectorXd Strain = Elements[eTag]->GetVTKResponse("Strain");
            OutputFile << Strain(0) << " " << Strain(1) << " " << Strain(2) << " " << Strain(3) << " " << Strain(4) << " " << Strain(5) << "\n";
        }

        //Writes the acceleration at each node.
        OutputFile << "\nStresses 6 " << Elements.size() << " float\n";
        for(auto it : Elements){
            auto &eTag = it.first;

            //Gets the material stress.
            Eigen::VectorXd Stress = Elements[eTag]->GetVTKResponse("Stress");
            OutputFile << Stress(0) << " " << Stress(1) << " " << Stress(2) << " "  << Stress(3) << " " << Stress(4) << " " << Stress(5) << "\n";
        }

        OutputFile.close();
        Counter++;
    }
}

//Fix blank spaces provided by user in path.
std::string 
Recorder::GetSpacedName(std::string theFile, std::string toReplace){
    //Auxiliar variables.
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

//Creates the directory to store results.
void 
Recorder::MakeFolder(std::string dirname){
    //Creates the folder where solutions are stored.
    std::string firstPath = GetSpacedName(filePath, "\\ ");
    std::string theFolder = firstPath + dirname + Combo;
    mkdir(theFolder.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}
