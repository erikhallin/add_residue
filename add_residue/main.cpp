#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <ctime>
#include "vec.h"

using namespace std;

float _deg2rad=0.0174532925;

struct st_pos
{
    st_pos()
    {
        x=y=z=0;
    }
    float x,y,z;
    string type;
};

vector<st_pos> vec_atoms_to_add;

st_pos get_pos_from_line(string line);
vec rotatePointAboutLine(vec p,float theta,vec p1,vec p2);
string force_3_digits(string input);
char from3to1(string letter_code);
bool load_rotamer(char residue_type);

int main( int argc, char *argv[] )
{
    cout<<"Software to add C-terminal residues of choice\n\n";

    bool have_seq_file=false;
    bool have_ipdb_file=false;
    bool have_opdb_file=false;
    bool have_residue_to_add=false;
    bool remove_pdb_input_h_and_o=true;

    string seq_file_name;
    string ipdb_file_name;
    string opdb_file_name;
    //string input_filename("input.pdb");//override with argument
    char residue_to_add='?';//override with argument

    for(int i=1;i<argc;i++)
    {
        //help
        if(string(argv[i])=="-help"||string(argv[i])=="-h"||string(argv[i])=="/help"||string(argv[i])=="/h")
        {
            /*cout<<"Instructions:\nFirst argument:  Input pdb structure file\nSecond argument: The residue to add in one-letter code\n\n";
            cout<<"Example: add_residue.exe my_structure.pdb A\nThis will add one alanine to your pdb file\n\n";
            cout<<"Give a filename that does not exist to make a new file with only that residue\n\n";*/

            cout<<"Instructions\n";
            cout<<"-h    :  Show help\n";
            cout<<"-seq  :  Name of file containing the protein sequence\n";
            cout<<"-ipdb :  Name of the input pdb structure file which to the residue should be added\n";
            cout<<"-opdb :  Name of the output pdb structure file having the residue added\n";
            cout<<"-res  :  Residue to add if not sequence file is used\n\n";

            return 9;
        }

        //sequence file
        if(string(argv[i])=="-seq")
        {
            //check arg after else error
            if(argc<=i-1 || argv[i+1][0]=='-')
            {
                cout<<"ERROR: Missing argument after -seq\n";
                return 9;
            }
            seq_file_name=argv[i+1];
            have_seq_file=true;
            i++;//skip reading next argument;
            continue;
        }

        //input pdb file
        if(string(argv[i])=="-ipdb")
        {
            //check arg after else error
            if(argc<=i-1 || argv[i+1][0]=='-')
            {
                cout<<"ERROR: Missing argument after -ipdb\n";
                return 9;
            }
            ipdb_file_name=argv[i+1];
            have_ipdb_file=true;
            i++;//skip reading next argument;
            continue;
        }

        //output pdb file
        if(string(argv[i])=="-opdb")
        {
            //check arg after else error
            if(argc<=i-1 || argv[i+1][0]=='-')
            {
                cout<<"ERROR: Missing argument after -opdb\n";
                return 9;
            }
            opdb_file_name=argv[i+1];
            have_opdb_file=true;
            i++;//skip reading next argument;
            continue;
        }

        //residue to add
        if(string(argv[i])=="-res")
        {
            //check arg after else error
            if(argc<=i-1 || argv[i+1][0]=='-')
            {
                cout<<"ERROR: Missing argument after -pdb\n";
                return 9;
            }
            string input_string(argv[i+1]);
            if((int)input_string.size()>1)
            {
                //convert from 3 letter to 1
                residue_to_add=from3to1(input_string);

                if(from3to1(input_string)=='?')
                {
                    cout<<"ERROR: Input residue is unknown: "<<input_string<<endl;
                    return 9;
                }
            }
            //otherwise take first char
            else residue_to_add=argv[i+1][0];

            //to upper case
            if(residue_to_add<=122 && residue_to_add>=97) residue_to_add-=32;

            have_residue_to_add=true;
            i++;//skip reading next argument;
            continue;
        }

        //do not remove hydrogens or carboxyl oxygen 2
        if(string(argv[i])=="-notremove")
        {
            remove_pdb_input_h_and_o=false;
        }

        //argument unstated, ignored
        cout<<"ERROR: Unknown argument: "<<argv[i]<<endl;
        return 9;


        /*switch(i)
        {
            //case 0://software name

            case 1://input filename
            {
                ipdb_file_name=argv[1];
            }break;

            case 2://input filename
            {
                residue_to_add=argv[2][0];

                //to upper case
                if(residue_to_add<=122 && residue_to_add>=97) residue_to_add-=32;
            }break;
        }*/
    }

    string line,word;
    srand(time(0));

    //status from input arguments
    //if no input pdb file, create empty as long as there is a seq file or res to add
    if(!have_ipdb_file)
    {
        //make new file...
        if(have_opdb_file) ipdb_file_name=opdb_file_name;
        else ipdb_file_name="tmp.pdb";

        if(!have_seq_file && !have_residue_to_add)
        {
            cout<<"ERROR: No sequence file or residue to add have been declared\n";
            return 10;
        }
    }
    //if no output file declared, make one up
    if(!have_opdb_file)
    {
        //add "_out" after file name
        string new_file_name(ipdb_file_name.c_str(),(int)ipdb_file_name.size()-4);
        new_file_name.append("_out.pdb");
        opdb_file_name=new_file_name;
    }
    //if input and output have same name, ???
    if(have_ipdb_file && have_opdb_file && ipdb_file_name==opdb_file_name)
    {
        //make backup of input file
        string backup_file_name(ipdb_file_name.c_str(),(int)ipdb_file_name.size()-4);
        backup_file_name.append("_backup.pdb");
        ofstream backup_file(backup_file_name.c_str());
        if(backup_file==0)
        {
            cout<<"ERROR: Could not create input backup file\n";
            return 11;
        }
        ifstream input_for_backup_file(ipdb_file_name.c_str());
        if(input_for_backup_file==0)
        {
            cout<<"ERROR: Could not find input pdb file\n";
            return 12;
        }
        while(getline(input_for_backup_file,line))
        {
            backup_file<<line<<endl;
        }
        input_for_backup_file.close();
        backup_file.close();

        //input file name have now been updated
        ipdb_file_name=backup_file_name;
    }
    //if no sequence or residue input, exit
    if(!have_seq_file && residue_to_add=='?')
    {
        cout<<"ERROR: No sequence file or specified residue to add have been declared\n";
        return 9;
    }

    //read seq file
    string input_sequence;
    if(have_seq_file)
    {
        ifstream seq_file(seq_file_name.c_str());
        if(seq_file==0)
        {
            cout<<"ERROR: Sequence file could not be found\n";
            return 7;
        }

        //char current_aa=pdb_seq[0];
        while(getline(seq_file,line))
        {
            if(line[0]=='>') continue; //skip header

            //save line
            input_sequence.append(line);
        }
    }

    //vector<st_pos> vec_atoms_to_add;

    //read pdb
    ifstream input_file(ipdb_file_name.c_str());
    if(input_file==0)
    {
        //create a new file then exit
        cout<<"Input pdb file not found, making a new one\n";
        cout<<"Creating new pdb file with one residue: "<<opdb_file_name<<endl;

        /*string output_filename;
        output_filename=string(ipdb_file_name.c_str(),(int)ipdb_file_name.size()-4);
        output_filename.append(1,'_');
        output_filename.append(1,residue_to_add);
        output_filename.append(".pdb");
        ofstream output_file(output_filename.c_str());*/

        ofstream output_file(opdb_file_name.c_str());
        if(output_file==0)
        {
            cout<<"ERROR: Could not create output file\n";
            return 5;
        }

        //have to know which residue to add
        if(have_seq_file) residue_to_add=input_sequence[0];

        //have to read rotamer library
        if(!load_rotamer(residue_to_add))
        {
            cout<<"ERROR: Could not load current rotamer\n";
            return 16;
        }

        int atom_counter=0;
        for(int i=0;i<(int)vec_atoms_to_add.size();i++)
        {
            stringstream ss_resnum;
            ss_resnum<<1;
            string resnum(ss_resnum.str());
            stringstream ss_atomnum;
            ss_atomnum<<++atom_counter;
            string atomnum(ss_atomnum.str());

            //pos
            stringstream ss_x;
            ss_x<<vec_atoms_to_add[i].x;
            string pos_x(ss_x.str());
            pos_x=force_3_digits(pos_x);

            stringstream ss_y;
            ss_y<<vec_atoms_to_add[i].y;
            string pos_y(ss_y.str());
            pos_y=force_3_digits(pos_y);

            stringstream ss_z;
            ss_z<<vec_atoms_to_add[i].z;
            string pos_z(ss_z.str());
            pos_z=force_3_digits(pos_z);

            string new_line("ATOM      ?  ??????? A   1       0.000   0.000   0.000  1.00  0.00            ");
            new_line[77]=vec_atoms_to_add[i].type[0];
            //update residue number
            for(int c=(int)resnum.size()-1;c>-1;c--)
            {
                new_line[27-(int)resnum.size()-1+c]=resnum[c];
            }

            //update atom number
            for(int c=(int)atomnum.size()-1;c>-1;c--)
            {
                new_line[12-(int)atomnum.size()-1+c]=atomnum[c];
            }

            //update type
            for(int c=(int)vec_atoms_to_add[i].type.size()-1;c>-1;c--)
            {
                new_line[21-(int)vec_atoms_to_add[i].type.size()-1+c]=vec_atoms_to_add[i].type[c];
            }

            //update pos
            for(int c=(int)pos_x.size()-1;c>-1;c--)
            {
                new_line[39-(int)pos_x.size()-1+c]=pos_x[c];
            }
            for(int c=(int)pos_y.size()-1;c>-1;c--)
            {
                new_line[47-(int)pos_y.size()-1+c]=pos_y[c];
            }
            for(int c=(int)pos_z.size()-1;c>-1;c--)
            {
                new_line[55-(int)pos_z.size()-1+c]=pos_z[c];
            }

            output_file<<new_line<<endl;
        }
        output_file.close();

        cout<<"Complete\n";
        return 0;


        //cout<<"ERROR: Could not find input file\n";
        //return 1;
    }

    //count residues (NB! only 999 residues possible)
    //and store all atom pos
    //and get pdb sequence
    vector<st_pos> vec_atoms_all;
    int end_residue_number=0;
    int end_atom_number=0;
    string pdb_seq;
    bool warning_odd_sequence=false;
    string last_line;
    int carboxyl_o_removed_counter=0;
    while(getline(input_file,line))
    {
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            //cout<<line<<endl;

            last_line=line;
            string number_as_text;
            number_as_text.append(1,line[23]);
            number_as_text.append(1,line[24]);
            number_as_text.append(1,line[25]);
            int res_number=(int)atof(number_as_text.c_str());
            if(res_number>end_residue_number)
            {
                end_residue_number=res_number;

                //save new residue to pdb seq file
                string current_aa3;
                current_aa3.append(1,line[17]);
                current_aa3.append(1,line[18]);
                current_aa3.append(1,line[19]);
                pdb_seq.append( 1, from3to1(current_aa3) );

                if(from3to1(current_aa3)=='?') warning_odd_sequence=true;
            }

            number_as_text="";
            number_as_text.append(1,line[5]);
            number_as_text.append(1,line[6]);
            number_as_text.append(1,line[7]);
            number_as_text.append(1,line[8]);
            number_as_text.append(1,line[9]);
            number_as_text.append(1,line[10]);
            int atom_number=(int)atof(number_as_text.c_str());
            if(atom_number>end_atom_number) end_atom_number=atom_number;

            //do not count carboxyl oxygen2
            if(remove_pdb_input_h_and_o&&line[13]=='O'&&line[14]=='2')
            {
                carboxyl_o_removed_counter++;
                end_atom_number--;
                continue;
            }

            //save atom pos and type
            vec_atoms_all.push_back(get_pos_from_line(line));
            vec_atoms_all.back().type=line[13];
        }
    }
    cout<<"Residues in protein: "<<end_residue_number<<endl;
    cout<<"Atoms in protein: "<<end_atom_number<<endl;
    if(remove_pdb_input_h_and_o) cout<<"C-terminal oxygens removed: "<<carboxyl_o_removed_counter<<endl;
    if(end_residue_number==0)
    {
        cout<<"ERROR: No residues found in input file\n";
        return 2;
    }
    if(warning_odd_sequence)
    {
        cout<<"WARNING: The sequence file contains unknown residue types\n";
    }


    /*//get seq from pdb
    string pdb_seq;
    string current_aa3("XXX");
    bool warning_odd_sequence=false;
    input_file.clear();
    input_file.seekg(0);
    while(getline(input_file,line))
    {
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            if(line[17]!=current_aa3[0]||line[18]!=current_aa3[1]||line[19]!=current_aa3[2])
            {
                current_aa3="";
                current_aa3.append(1,line[17]);
                current_aa3.append(1,line[18]);
                current_aa3.append(1,line[19]);

                pdb_seq.append( 1, from3to1(current_aa3) );

                if(from3to1(current_aa3)=='?') warning_odd_sequence=true;
            }
        }
    }
    if(warning_odd_sequence)
    {
        cout<<"WARNING: The sequence file contains unknown residue types\n";
    }*/



    if(!pdb_seq.empty() && have_seq_file)
    {
        //find next residue to add
        int last_matching_aa_ind=-1;
        for(int i=0;i<(int)input_sequence.size() && i<(int)pdb_seq.size();i++)
        {
            if(input_sequence[i]==pdb_seq[i]) last_matching_aa_ind=i;
        }
        //check if mismatch in sequence
        if(last_matching_aa_ind==-1 || last_matching_aa_ind!=(int)pdb_seq.size()-1)
        {
            cout<<"ERROR: Sequence file does not match with the structure sequence\n";
            return 8;
        }
        //check if end of input sequence have been reached
        if(last_matching_aa_ind==(int)input_sequence.size()-1)
        {
            cout<<"ERROR: The end of the input sequence have been reached, whole structure is already build\n";
            return 9;
        }

        residue_to_add=input_sequence[last_matching_aa_ind+1];
    }
    else
    {
        //start new pdb file
        if(have_seq_file) residue_to_add=input_sequence[0];

        //residue type from argument instead
        if(residue_to_add=='?')
        {
            cout<<"ERROR: Residue to add have not been specified\n";
            return 9;
        }
    }


    cout<<"Input pdb file: "<<ipdb_file_name<<endl;
    cout<<"Residue to be added: "<<residue_to_add<<endl;


    //load rotamer
    if(!load_rotamer(residue_to_add))
    {
        cout<<"ERROR: Could not load current rotamer\n";
        return 17;
    }


    //find c-term
    input_file.clear();
    input_file.seekg(0);
    st_pos pos_prev_n;
    st_pos pos_prev_ca;
    st_pos pos_prev_c;
    st_pos pos_prev_o;
    while(getline(input_file,line))
    {
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            //cout<<line<<endl;

            //find end residue
            string number_as_text;
            number_as_text.append(1,line[23]);
            number_as_text.append(1,line[24]);
            number_as_text.append(1,line[25]);
            int res_number=(int)atof(number_as_text.c_str());
            if(res_number==end_residue_number)
            {
                //save backbone pos
                if(line[13]=='N'&&line[14]==' ')
                {
                    pos_prev_n=get_pos_from_line(line);
                }
                if(line[13]=='C'&&line[14]=='A')
                {
                    pos_prev_ca=get_pos_from_line(line);
                }
                if(line[13]=='C'&&line[14]==' ')
                {
                    pos_prev_c=get_pos_from_line(line);
                }
                if((line[13]=='O'&&line[14]==' ') || (line[13]=='O'&&line[14]=='1'))
                {
                    pos_prev_o=get_pos_from_line(line);
                }
            }
        }
    }

    //calc new residue pos
    st_pos prev_n_to_ca_dir;
    prev_n_to_ca_dir.x=pos_prev_ca.x-pos_prev_n.x;
    prev_n_to_ca_dir.y=pos_prev_ca.y-pos_prev_n.y;
    prev_n_to_ca_dir.z=pos_prev_ca.z-pos_prev_n.z;
    //normalize
    vec dir(prev_n_to_ca_dir.x,prev_n_to_ca_dir.y,prev_n_to_ca_dir.z);
    float length=dir.length();
    float bond_length=1.33;
    prev_n_to_ca_dir.x=prev_n_to_ca_dir.x/length;
    prev_n_to_ca_dir.y=prev_n_to_ca_dir.y/length;
    prev_n_to_ca_dir.z=prev_n_to_ca_dir.z/length;
    st_pos new_res_pos_n;
    new_res_pos_n.x=pos_prev_c.x+prev_n_to_ca_dir.x*bond_length;
    new_res_pos_n.y=pos_prev_c.y+prev_n_to_ca_dir.y*bond_length;
    new_res_pos_n.z=pos_prev_c.z+prev_n_to_ca_dir.z*bond_length;
    //cout<<pos_prev_n.x<<pos_prev_n.y<<pos_prev_n.z<<endl;
    //cout<<pos_prev_ca.x<<pos_prev_ca.y<<pos_prev_ca.z<<endl;
    //cout<<pos_prev_c.x<<pos_prev_c.y<<pos_prev_c.z<<endl;
    //cout<<prev_n_to_ca_dir.x<<prev_n_to_ca_dir.y<<prev_n_to_ca_dir.z<<endl;

    //guess rotation of new residue to find pos furthest away from prev aa
    vector<vec> vec_rotations;
    vector<int> vec_clashcount;
    vector<float> vec_distance;
    int total_tries=1000;
    float longest_dist=0;
    vec guessed_rot_longest_dist;
    vec guessed_rot_longest_dist_lowest_clash;
    for(int i=0;i<total_tries;i++)
    {
        vec start_pos(0.551, -1.198, -0.766);
        //rotate
        float rot_x=rand()%360;
        float rot_y=rand()%360;
        float rot_z=rand()%360;
        //guessed_rot=vec(rot_x,rot_y,rot_z);
        vec_rotations.push_back(vec(rot_x,rot_y,rot_z));

        vec rotated_pos=start_pos;
        rotated_pos=rotatePointAboutLine(rotated_pos,rot_x,vec(0,0,0),vec(1,0,0));
        rotated_pos=rotatePointAboutLine(rotated_pos,rot_y,vec(0,0,0),vec(0,1,0));
        rotated_pos=rotatePointAboutLine(rotated_pos,rot_z,vec(0,0,0),vec(0,0,1));

        //add absolute pos
        rotated_pos.x+=new_res_pos_n.x;
        rotated_pos.y+=new_res_pos_n.y;
        rotated_pos.z+=new_res_pos_n.z;

        //measure distance
        vec prev_c_dir(rotated_pos.x-pos_prev_c.x, rotated_pos.y-pos_prev_c.y, rotated_pos.z-pos_prev_c.z);
        float dist=prev_c_dir.length();
        vec_distance.push_back(dist);

        if(dist>longest_dist)
        {
            longest_dist=dist;
            guessed_rot_longest_dist=vec(rot_x,rot_y,rot_z);
        }

        //calc final pos of all atoms for the new residue for the clash test
        vector<vec> vec_pos_residue;
        for(int ind=0;ind<(int)vec_atoms_to_add.size();ind++)
        {
            //rotate
            vec new_pos(vec_atoms_to_add[ind].x,vec_atoms_to_add[ind].y,vec_atoms_to_add[ind].z);
            new_pos=rotatePointAboutLine(new_pos,vec_rotations.back().x,vec(0,0,0),vec(1,0,0));
            new_pos=rotatePointAboutLine(new_pos,vec_rotations.back().y,vec(0,0,0),vec(0,1,0));
            new_pos=rotatePointAboutLine(new_pos,vec_rotations.back().z,vec(0,0,0),vec(0,0,1));
            //translate
            new_pos.x+=new_res_pos_n.x;
            new_pos.y+=new_res_pos_n.y;
            new_pos.z+=new_res_pos_n.z;

            vec_pos_residue.push_back(new_pos);
        }

        //clash test
        float clash_dist_threshold=1.0;
        float clash_trig2=clash_dist_threshold*clash_dist_threshold;
        int clash_counter=0;
        for(int i_all=0;i_all<(int)vec_atoms_all.size();i_all++)
        {
            for(int i_res=0;i_res<(int)vec_pos_residue.size();i_res++)
            {
                float dist2=(vec_atoms_all[i_all].x-vec_pos_residue[i_res].x)*(vec_atoms_all[i_all].x-vec_pos_residue[i_res].x)+
                            (vec_atoms_all[i_all].y-vec_pos_residue[i_res].y)*(vec_atoms_all[i_all].y-vec_pos_residue[i_res].y)+
                            (vec_atoms_all[i_all].z-vec_pos_residue[i_res].z)*(vec_atoms_all[i_all].z-vec_pos_residue[i_res].z);
                if(dist2<clash_trig2) clash_counter++;
            }
        }
        vec_clashcount.push_back(clash_counter);
        //cout<<clash_counter<<endl;
    }

    //select suitable rotation alternative
    int clash_tolerance=0;
    while(true)
    {
        vector<int> vec_set_w_current_clash_count;
        for(int i=0;i<(int)vec_clashcount.size();i++)
        {
            if(vec_clashcount[i]<=clash_tolerance) vec_set_w_current_clash_count.push_back(i);
        }

        if(vec_set_w_current_clash_count.empty())
        {
            clash_tolerance++;
        }
        else
        {
            //select alternative with longest distance
            float longest_dist_in_set=0;
            int longest_ind=0;
            for(int i=0;i<(int)vec_set_w_current_clash_count.size();i++)
            {
                if(vec_distance[ vec_set_w_current_clash_count[i] ]>longest_dist_in_set)
                {
                    longest_dist_in_set=vec_distance[ vec_set_w_current_clash_count[i] ];
                    longest_ind=vec_set_w_current_clash_count[i];
                }
            }

            //save best choice
            guessed_rot_longest_dist_lowest_clash=vec_rotations[longest_ind];

            break;
        }
    }

    //move res into place
    for(int i=0;i<(int)vec_atoms_to_add.size();i++)
    {
        //rotate
        vec new_pos(vec_atoms_to_add[i].x,vec_atoms_to_add[i].y,vec_atoms_to_add[i].z);
        new_pos=rotatePointAboutLine(new_pos,guessed_rot_longest_dist_lowest_clash.x,vec(0,0,0),vec(1,0,0));
        new_pos=rotatePointAboutLine(new_pos,guessed_rot_longest_dist_lowest_clash.y,vec(0,0,0),vec(0,1,0));
        new_pos=rotatePointAboutLine(new_pos,guessed_rot_longest_dist_lowest_clash.z,vec(0,0,0),vec(0,0,1));
        //translate
        new_pos.x+=new_res_pos_n.x;
        new_pos.y+=new_res_pos_n.y;
        new_pos.z+=new_res_pos_n.z;
        //cap near zero
        if(fabs(new_pos.x)<0.001) new_pos.x=0;
        if(fabs(new_pos.y)<0.001) new_pos.y=0;
        if(fabs(new_pos.z)<0.001) new_pos.z=0;
        //cap 3 digits
        new_pos.x=float(int((new_pos.x+0.0005)*1000.0))/1000.0;
        new_pos.y=float(int((new_pos.y+0.0005)*1000.0))/1000.0;
        new_pos.z=float(int((new_pos.z+0.0005)*1000.0))/1000.0;

        vec_atoms_to_add[i].x=new_pos.x;
        vec_atoms_to_add[i].y=new_pos.y;
        vec_atoms_to_add[i].z=new_pos.z;
    }

    //build residue in new file
    /*string output_filename("output.pdb");
    if(ipdb_file_name!="input.pdb")
    {
        output_filename=string(ipdb_file_name.c_str(),(int)ipdb_file_name.size()-4);
        output_filename.append(1,'_');
        output_filename.append(1,residue_to_add);
        output_filename.append(".pdb");
    }*/
    cout<<"Creating new pdb file with added residue: "<<opdb_file_name<<endl;
    ofstream output_file(opdb_file_name.c_str());
    if(output_file==0)
    {
        cout<<"ERROR: Could not create output file\n";
        return 4;
    }
    input_file.clear();
    input_file.seekg(0);
    bool final_pos_found=false;
    bool new_residue_added=false;
    char chain_name=' ';//use same as for previous residue
    bool add_atom_type=false;
    while(getline(input_file,line))
    {
        if(!final_pos_found)
        {
            if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
            {
                //find end residue
                string number_as_text;
                number_as_text.append(1,line[23]);
                number_as_text.append(1,line[24]);
                number_as_text.append(1,line[25]);
                int res_number=(int)atof(number_as_text.c_str());
                if(res_number==end_residue_number)
                {
                    //change name of C-term oxygen if O1, should be only O
                    if(line[13]=='O'&&line[14]=='1') line[14]=' ';
                }
                //find end atom
                number_as_text="";
                number_as_text.append(1,line[5]);
                number_as_text.append(1,line[6]);
                number_as_text.append(1,line[7]);
                number_as_text.append(1,line[8]);
                number_as_text.append(1,line[9]);
                number_as_text.append(1,line[10]);
                int atom_number=(int)atof(number_as_text.c_str());
                //cout<<atom_number<<" : "<<end_atom_number<<endl;
                if(atom_number==end_atom_number)
                {
                    //cout<<"Final atom number: "<<atom_number<<endl;

                    final_pos_found=true;
                    //store chain name
                    chain_name=line[21];
                    //check if atom type is stated
                    if((int)line.size()>=78)
                    {
                        if(line[77]!=' ') add_atom_type=true;
                    }

                    //save last line
                    output_file<<line<<endl;

                    //cout<<line<<endl;
                }
            }
        }
        //if extra text in the file, add text now
        else if(!new_residue_added/*&&(line[0]!='A'||line[1]!='T'||line[2]!='O'||line[3]!='M')*/)
        {
            //add new residue before extra text from original file
            new_residue_added=true;
            //add_new_residue_lines(vec_atoms_to_add,output_file,end_residue_number,end_atom_number);

            for(int i=0;i<(int)vec_atoms_to_add.size();i++)
            {
                stringstream ss_resnum;
                ss_resnum<<end_residue_number+1;
                string resnum(ss_resnum.str());
                stringstream ss_atomnum;
                ss_atomnum<<++end_atom_number;
                string atomnum(ss_atomnum.str());

                //pos
                stringstream ss_x;
                //ss_x.precision(3);
                ss_x<<vec_atoms_to_add[i].x;
                //ss_x<<int(vec_atoms_to_add[i].x)<<'.'<<fabs(int(vec_atoms_to_add[i].x*1000.0)-int(vec_atoms_to_add[i].x)*1000.0);
                string pos_x(ss_x.str());
                pos_x=force_3_digits(pos_x);
                /*while((int)pos_x.size()<4) pos_x.append(1,'0');

                pos_x.append(1,'0');
                pos_x[(int)pos_x.size()-1]=pos_x[(int)pos_x.size()-2];
                pos_x[(int)pos_x.size()-2]=pos_x[(int)pos_x.size()-3];
                pos_x[(int)pos_x.size()-3]=pos_x[(int)pos_x.size()-4];
                pos_x[(int)pos_x.size()-4]='.';*/

                stringstream ss_y;
                //ss_y.precision(3);
                ss_y<<vec_atoms_to_add[i].y;
                //ss_y<<int(vec_atoms_to_add[i].y)<<'.'<<fabs(int(vec_atoms_to_add[i].y*1000.0)-int(vec_atoms_to_add[i].y)*1000.0);
                string pos_y(ss_y.str());
                pos_y=force_3_digits(pos_y);
                /*while((int)pos_y.size()<4) pos_y.append(1,'0');
                cout<<vec_atoms_to_add[i].y<<endl;
                cout<<int(vec_atoms_to_add[i].y*1000.0)<<endl;
                cout<<pos_y<<endl;
                pos_y.append(1,'0');
                cout<<pos_y<<endl;
                pos_y[(int)pos_y.size()-1]=pos_y[(int)pos_y.size()-2];
                pos_y[(int)pos_y.size()-2]=pos_y[(int)pos_y.size()-3];
                pos_y[(int)pos_y.size()-3]=pos_y[(int)pos_y.size()-4];
                pos_y[(int)pos_y.size()-4]='.';
                cout<<pos_y<<endl;*/

                stringstream ss_z;
                //ss_z.precision(3);
                ss_z<<vec_atoms_to_add[i].z;
                //ss_z<<int(vec_atoms_to_add[i].z)<<'.'<<fabs(int(vec_atoms_to_add[i].z*1000.0)-int(vec_atoms_to_add[i].z)*1000.0);
                string pos_z(ss_z.str());
                pos_z=force_3_digits(pos_z);
                /*while((int)pos_z.size()<4) pos_z.append(1,'0');
                pos_z.append(1,'0');
                pos_z[(int)pos_z.size()-1]=pos_z[(int)pos_z.size()-2];
                pos_z[(int)pos_z.size()-2]=pos_z[(int)pos_z.size()-3];
                pos_z[(int)pos_z.size()-3]=pos_z[(int)pos_z.size()-4];
                pos_z[(int)pos_z.size()-4]='.';*/

                //chain name pos marked with (A)
                string new_line("ATOM      ?  ??????? A   1       0.000   0.000   0.000  1.00  0.00            ");
                new_line[21]=chain_name;
                if(add_atom_type) new_line[77]=vec_atoms_to_add[i].type[0];

                //update residue number
                for(int c=(int)resnum.size()-1;c>-1;c--)
                {
                    new_line[27-(int)resnum.size()-1+c]=resnum[c];
                }

                //update atom number
                for(int c=(int)atomnum.size()-1;c>-1;c--)
                {
                    new_line[12-(int)atomnum.size()-1+c]=atomnum[c];
                }

                //update type
                for(int c=(int)vec_atoms_to_add[i].type.size()-1;c>-1;c--)
                {
                    new_line[21-(int)vec_atoms_to_add[i].type.size()-1+c]=vec_atoms_to_add[i].type[c];
                }

                //update pos
                for(int c=(int)pos_x.size()-1;c>-1;c--)
                {
                    new_line[39-(int)pos_x.size()-1+c]=pos_x[c];
                }
                for(int c=(int)pos_y.size()-1;c>-1;c--)
                {
                    new_line[47-(int)pos_y.size()-1+c]=pos_y[c];
                }
                for(int c=(int)pos_z.size()-1;c>-1;c--)
                {
                    new_line[55-(int)pos_z.size()-1+c]=pos_z[c];
                }

                //add to new file
                output_file<<new_line<<endl;

                //cout<<"orig: "<<line<<endl;
                //cout<<"newl: "<<new_line<<endl;
            }
        }

        //add line to new file
        //cout<<line<<endl;
        bool donotprint=false;
        if(final_pos_found && (line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')) donotprint=true;//ignore atoms after new residue added
        if(final_pos_found && !new_residue_added) donotprint=true;
        if(!donotprint) output_file<<line<<endl;
    }
    if(!new_residue_added)
    {
        //end of file reached, add extra residue to new file now
        //add_new_residue_lines(vec_atoms_to_add,output_file,end_residue_number,end_atom_number);

        for(int i=0;i<(int)vec_atoms_to_add.size();i++)
        {
            stringstream ss_resnum;
            ss_resnum<<end_residue_number+1;
            string resnum(ss_resnum.str());
            stringstream ss_atomnum;
            ss_atomnum<<++end_atom_number;
            string atomnum(ss_atomnum.str());

            //pos
            stringstream ss_x;
            //ss_x.precision(3);
            ss_x<<vec_atoms_to_add[i].x;
            //ss_x<<int(vec_atoms_to_add[i].x)<<'.'<<fabs(int(vec_atoms_to_add[i].x*1000.0)-int(vec_atoms_to_add[i].x)*1000.0);
            string pos_x(ss_x.str());
            pos_x=force_3_digits(pos_x);
            //cout<<"this: "<<vec_atoms_to_add[i].x<<endl;
            /*while((int)pos_x.size()<4) pos_x.append(1,'0');

            pos_x.append(1,'0');
            pos_x[(int)pos_x.size()-1]=pos_x[(int)pos_x.size()-2];
            pos_x[(int)pos_x.size()-2]=pos_x[(int)pos_x.size()-3];
            pos_x[(int)pos_x.size()-3]=pos_x[(int)pos_x.size()-4];
            pos_x[(int)pos_x.size()-4]='.';*/

            stringstream ss_y;
            //ss_y.precision(3);
            ss_y<<vec_atoms_to_add[i].y;
            //ss_y<<int(vec_atoms_to_add[i].y)<<'.'<<fabs(int(vec_atoms_to_add[i].y*1000.0)-int(vec_atoms_to_add[i].y)*1000.0);
            string pos_y(ss_y.str());
            pos_y=force_3_digits(pos_y);
            /*while((int)pos_y.size()<4) pos_y.append(1,'0');
            cout<<vec_atoms_to_add[i].y<<endl;
            cout<<int(vec_atoms_to_add[i].y*1000.0)<<endl;
            cout<<pos_y<<endl;
            pos_y.append(1,'0');
            cout<<pos_y<<endl;
            pos_y[(int)pos_y.size()-1]=pos_y[(int)pos_y.size()-2];
            pos_y[(int)pos_y.size()-2]=pos_y[(int)pos_y.size()-3];
            pos_y[(int)pos_y.size()-3]=pos_y[(int)pos_y.size()-4];
            pos_y[(int)pos_y.size()-4]='.';
            cout<<pos_y<<endl;*/

            stringstream ss_z;
            //ss_z.precision(3);
            ss_z<<vec_atoms_to_add[i].z;
            //ss_z<<int(vec_atoms_to_add[i].z)<<'.'<<fabs(int(vec_atoms_to_add[i].z*1000.0)-int(vec_atoms_to_add[i].z)*1000.0);
            string pos_z(ss_z.str());
            pos_z=force_3_digits(pos_z);

            //chain name pos marked with (A)
            string new_line("ATOM      ?  ??????? A   1       0.000   0.000   0.000  1.00  0.00            ");
            new_line[21]=chain_name;
            if(add_atom_type) new_line[77]=vec_atoms_to_add[i].type[0];

            //update residue number
            for(int c=(int)resnum.size()-1;c>-1;c--)
            {
                new_line[27-(int)resnum.size()-1+c]=resnum[c];
            }

            //update atom number
            for(int c=(int)atomnum.size()-1;c>-1;c--)
            {
                new_line[12-(int)atomnum.size()-1+c]=atomnum[c];
            }

            //update type
            for(int c=(int)vec_atoms_to_add[i].type.size()-1;c>-1;c--)
            {
                new_line[21-(int)vec_atoms_to_add[i].type.size()-1+c]=vec_atoms_to_add[i].type[c];
            }

            //update pos
            for(int c=(int)pos_x.size()-1;c>-1;c--)
            {
                new_line[39-(int)pos_x.size()-1+c]=pos_x[c];
            }
            for(int c=(int)pos_y.size()-1;c>-1;c--)
            {
                new_line[47-(int)pos_y.size()-1+c]=pos_y[c];
            }
            for(int c=(int)pos_z.size()-1;c>-1;c--)
            {
                new_line[55-(int)pos_z.size()-1+c]=pos_z[c];
            }

            //add to new file
            output_file<<new_line<<endl;

            //cout<<"orig: "<<line<<endl;
            //cout<<"newl: "<<new_line<<endl;
        }
    }

    input_file.close();
    output_file.close();

    //add c-term oxygen?

    cout<<"Complete\n";

    return 0;
}

st_pos get_pos_from_line(string line)
{
    st_pos pos;
    string number_as_text;
    number_as_text.append(1,line[30]);
    number_as_text.append(1,line[31]);
    number_as_text.append(1,line[32]);
    number_as_text.append(1,line[33]);
    number_as_text.append(1,line[34]);
    number_as_text.append(1,line[35]);
    number_as_text.append(1,line[36]);
    number_as_text.append(1,line[37]);
    pos.x=atof(number_as_text.c_str());
    number_as_text="";
    number_as_text.append(1,line[38]);
    number_as_text.append(1,line[39]);
    number_as_text.append(1,line[40]);
    number_as_text.append(1,line[41]);
    number_as_text.append(1,line[42]);
    number_as_text.append(1,line[43]);
    number_as_text.append(1,line[44]);
    number_as_text.append(1,line[45]);
    pos.y=atof(number_as_text.c_str());
    number_as_text="";
    number_as_text.append(1,line[46]);
    number_as_text.append(1,line[47]);
    number_as_text.append(1,line[48]);
    number_as_text.append(1,line[49]);
    number_as_text.append(1,line[50]);
    number_as_text.append(1,line[51]);
    number_as_text.append(1,line[52]);
    number_as_text.append(1,line[53]);
    pos.z=atof(number_as_text.c_str());

    return pos;
}

vec rotatePointAboutLine(vec p,float theta,vec p1,vec p2)
{
   //cout<<"p: "<<p.x<<", "<<p.y<<", "<<p.z<<" Theta: "<<theta<<endl;

   //Rotation of a point in 3 dimensional space by theta about an arbitrary axes defined by a line
   //between two points P1 = (x1,y1,z1) and P2 = (x2,y2,z2) can be achieved by the following steps
   vec u,q1,q2;
   float d;

   // Step 1 translate space so that the rotation axis passes through the origin
   q1.x = p.x - p1.x;
   q1.y = p.y - p1.y;
   q1.z = p.z - p1.z;

   u.x = p2.x - p1.x;
   u.y = p2.y - p1.y;
   u.z = p2.z - p1.z;
   //Normalise(&u);
   u=u.unit();
   d = sqrt(u.y*u.y + u.z*u.z);

   // Step 2 rotate space about the x axis so that the rotation axis lies in the xz plane
   if (d != 0)
   {
      q2.x = q1.x;
      q2.y = q1.y * u.z / d - q1.z * u.y / d;
      q2.z = q1.y * u.y / d + q1.z * u.z / d;
   }
   else
   {
      q2 = q1;
   }

   // Step 3 rotate space about the y axis so that the rotation axis lies along the z axis
   q1.x = q2.x * d - q2.z * u.x;
   q1.y = q2.y;
   q1.z = q2.x * u.x + q2.z * d;

   // Step 4 perform the desired rotation by theta about the z axis
   q2.x = q1.x * cosf(theta*_deg2rad) - q1.y * sinf(theta*_deg2rad);
   q2.y = q1.x * sinf(theta*_deg2rad) + q1.y * cosf(theta*_deg2rad);
   q2.z = q1.z;

   // Inverse of step 3
   q1.x =   q2.x * d + q2.z * u.x;
   q1.y =   q2.y;
   q1.z = - q2.x * u.x + q2.z * d;

   // Inverse of step 2
   if (d != 0)
   {
      q2.x =   q1.x;
      q2.y =   q1.y * u.z / d + q1.z * u.y / d;
      q2.z = - q1.y * u.y / d + q1.z * u.z / d;
   }
   else
   {
      q2 = q1;
   }

   // Inverse of step 1
   q1.x = q2.x + p1.x;
   q1.y = q2.y + p1.y;
   q1.z = q2.z + p1.z;

   //cout<<"q1: "<<q1.x<<", "<<q1.y<<", "<<q1.z<<endl;

   return q1;
}

string force_3_digits(string input)
{
    int digit_counter=0;
    bool dot_found=false;
    for(int i=0;i<(int)input.size();i++)
    {
        if(dot_found) digit_counter++;

        if(!dot_found && input[i]=='.') dot_found=true;
    }
    //add dot if not found
    if(!dot_found) input.append(1,'.');
    //add 0 or remove digits
    if(digit_counter<3)
    {
        for(int i=digit_counter;i<3;i++)
        {
            input.append(1,'0');
        }
    }
    else if(digit_counter>3)
    {
        //cut end
        input=string(input.c_str(),(int)input.size()-digit_counter+3);
    }

    return input;
}

char from3to1(string letter_code)
{
	if(letter_code=="ALA"||letter_code=="Ala"||letter_code=="ala") return 'A';
	if(letter_code=="ARG"||letter_code=="Arg"||letter_code=="arg") return 'R';
	if(letter_code=="ASN"||letter_code=="Asn"||letter_code=="asn") return 'N';
	if(letter_code=="ASP"||letter_code=="Asp"||letter_code=="asp") return 'D';
	if(letter_code=="CYS"||letter_code=="Cys"||letter_code=="cys") return 'C';
	if(letter_code=="GLU"||letter_code=="Glu"||letter_code=="glu") return 'E';
	if(letter_code=="GLN"||letter_code=="Gln"||letter_code=="gln") return 'Q';
	if(letter_code=="GLY"||letter_code=="Gly"||letter_code=="gly") return 'G';
	if(letter_code=="HIS"||letter_code=="His"||letter_code=="his") return 'H';
	if(letter_code=="ILE"||letter_code=="Ile"||letter_code=="ile") return 'I';
	if(letter_code=="LEU"||letter_code=="Leu"||letter_code=="leu") return 'L';
	if(letter_code=="LYS"||letter_code=="Lys"||letter_code=="lys") return 'K';
	if(letter_code=="MET"||letter_code=="Met"||letter_code=="met") return 'M';
	if(letter_code=="PHE"||letter_code=="Phe"||letter_code=="phe") return 'F';
	if(letter_code=="PRO"||letter_code=="Pro"||letter_code=="pro") return 'P';
	if(letter_code=="SER"||letter_code=="Ser"||letter_code=="ser") return 'S';
	if(letter_code=="THR"||letter_code=="Thr"||letter_code=="thr") return 'T';
	if(letter_code=="TRP"||letter_code=="Trp"||letter_code=="trp") return 'W';
	if(letter_code=="TYR"||letter_code=="Tyr"||letter_code=="tyr") return 'Y';
	if(letter_code=="VAL"||letter_code=="Val"||letter_code=="val") return 'V';

	//unknown
	return '?';
}

bool load_rotamer(char residue_type)
{
    //read rotamer library
    cout<<"Loading rotamer library\n";
    string file_to_load_rotamer;
    switch(residue_type)
    {
        case 'A': file_to_load_rotamer="rotamers\\ala.pdb"; break;
        case 'R': file_to_load_rotamer="rotamers\\arg.pdb"; break;
        case 'N': file_to_load_rotamer="rotamers\\asn.pdb"; break;
        case 'D': file_to_load_rotamer="rotamers\\asp.pdb"; break;
        case 'C': file_to_load_rotamer="rotamers\\cys.pdb"; break;
        case 'Q': file_to_load_rotamer="rotamers\\gln.pdb"; break;
        case 'E': file_to_load_rotamer="rotamers\\glu.pdb"; break;
        case 'G': file_to_load_rotamer="rotamers\\gly.pdb"; break;
        case 'H': file_to_load_rotamer="rotamers\\his.pdb"; break;
        case 'I': file_to_load_rotamer="rotamers\\ile.pdb"; break;
        case 'L': file_to_load_rotamer="rotamers\\leu.pdb"; break;
        case 'K': file_to_load_rotamer="rotamers\\lys.pdb"; break;
        case 'M': file_to_load_rotamer="rotamers\\met.pdb"; break;
        case 'F': file_to_load_rotamer="rotamers\\phe.pdb"; break;
        case 'P': file_to_load_rotamer="rotamers\\pro.pdb"; break;
        case 'S': file_to_load_rotamer="rotamers\\ser.pdb"; break;
        case 'T': file_to_load_rotamer="rotamers\\thr.pdb"; break;
        case 'W': file_to_load_rotamer="rotamers\\trp.pdb"; break;
        case 'Y': file_to_load_rotamer="rotamers\\tyr.pdb"; break;
        case 'V': file_to_load_rotamer="rotamers\\val.pdb"; break;

        /*//if unix
        case 'A': file_to_load_rotamer="rotamers/ala.pdb"; break;
        case 'R': file_to_load_rotamer="rotamers/arg.pdb"; break;
        case 'N': file_to_load_rotamer="rotamers/asn.pdb"; break;
        case 'D': file_to_load_rotamer="rotamers/asp.pdb"; break;
        case 'C': file_to_load_rotamer="rotamers/cys.pdb"; break;
        case 'Q': file_to_load_rotamer="rotamers/gln.pdb"; break;
        case 'E': file_to_load_rotamer="rotamers/glu.pdb"; break;
        case 'G': file_to_load_rotamer="rotamers/gly.pdb"; break;
        case 'H': file_to_load_rotamer="rotamers/his.pdb"; break;
        case 'I': file_to_load_rotamer="rotamers/ile.pdb"; break;
        case 'L': file_to_load_rotamer="rotamers/leu.pdb"; break;
        case 'K': file_to_load_rotamer="rotamers/lys.pdb"; break;
        case 'M': file_to_load_rotamer="rotamers/met.pdb"; break;
        case 'F': file_to_load_rotamer="rotamers/phe.pdb"; break;
        case 'P': file_to_load_rotamer="rotamers/pro.pdb"; break;
        case 'S': file_to_load_rotamer="rotamers/ser.pdb"; break;
        case 'T': file_to_load_rotamer="rotamers/thr.pdb"; break;
        case 'W': file_to_load_rotamer="rotamers/trp.pdb"; break;
        case 'Y': file_to_load_rotamer="rotamers/tyr.pdb"; break;
        case 'V': file_to_load_rotamer="rotamers/val.pdb"; break;*/
    }
    ifstream file(file_to_load_rotamer.c_str());
    if(file==0)
    {
        cout<<"ERROR: Missing residue library file\n";
        return false;
    }
    //count rotamer alternatives
    int rotamer_total=0;
    string line;
    while(getline(file,line))
    {
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            string number_as_text;
            number_as_text.append(1,line[22]);
            number_as_text.append(1,line[23]);
            number_as_text.append(1,line[24]);
            number_as_text.append(1,line[25]);
            int rotamer_number=(int)atof(number_as_text.c_str());
            if(rotamer_number>rotamer_total) rotamer_total=rotamer_number;
        }
    }
    //select random rotamer
    int selected_rotamer_id=rand()%rotamer_total+1;
    cout<<"Selected rotamer: "<<selected_rotamer_id<<endl;
    //bool rotamer_found=false;
    file.clear();
    file.seekg(0);
    while(getline(file,line))
    {
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            //go to selected rotamer
            string number_as_text;
            number_as_text.append(1,line[22]);
            number_as_text.append(1,line[23]);
            number_as_text.append(1,line[24]);
            number_as_text.append(1,line[25]);
            int rotamer_number=(int)atof(number_as_text.c_str());
            if(rotamer_number!=selected_rotamer_id) continue;//skip this atom
            //rotamer_found=true;

            //if(line[25]=='2') break;//done with rotamer 1

            //ignore protons
            if(line[13]=='H') continue;
            //get pos
            vec_atoms_to_add.push_back(get_pos_from_line(line));
            //add type
            vec_atoms_to_add.back().type.append(1,line[13]);
            vec_atoms_to_add.back().type.append(1,line[14]);
            vec_atoms_to_add.back().type.append(1,line[15]);
            vec_atoms_to_add.back().type.append(1,line[16]);
            vec_atoms_to_add.back().type.append(1,line[17]);
            vec_atoms_to_add.back().type.append(1,line[18]);
            vec_atoms_to_add.back().type.append(1,line[19]);

            //center nitrogen
            vec_atoms_to_add.back().x+=1.458;
        }
    }

    return true;
}
