/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include <iostream>
#include <cerrno>

#define NO_ALL 1092


#ifdef WIN32
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
typedef unsigned short int uint16;
typedef int int32;
typedef short int int16;
typedef unsigned  int uint32;
typedef long long int64;
typedef unsigned  long long uint64;
typedef unsigned char uchar;
#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell
#define _TCHAR	char
#define _tmain	main

typedef unsigned short int uint16;
typedef int int32;
typedef short int int16;
typedef unsigned  int uint32;
typedef long long int64;
typedef unsigned long  long uint64;
typedef unsigned char uchar;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#endif

using namespace std;

char all_ind[NO_ALL][8] ={"HG00096", "HG00097", "HG00099", "HG00100", "HG00101", "HG00102", "HG00103", "HG00104", "HG00106", "HG00108", "HG00109", "HG00110", "HG00111", "HG00112", "HG00113", "HG00114", "HG00116", "HG00117", "HG00118", "HG00119", "HG00120", "HG00121", "HG00122", "HG00123", "HG00124", "HG00125", "HG00126", "HG00127", "HG00128", "HG00129", "HG00130", "HG00131", "HG00133", "HG00134", "HG00135", "HG00136", "HG00137", "HG00138", "HG00139", "HG00140", "HG00141", "HG00142", "HG00143", "HG00146", "HG00148", "HG00149", "HG00150", "HG00151", "HG00152", "HG00154", "HG00155", "HG00156", "HG00158", "HG00159", "HG00160", "HG00171", "HG00173", "HG00174", "HG00176", "HG00177", "HG00178", "HG00179", "HG00180", "HG00182", "HG00183", "HG00185", "HG00186", "HG00187", "HG00188", "HG00189", "HG00190", "HG00231", "HG00232", "HG00233", "HG00234", "HG00235", "HG00236", "HG00237", "HG00238", "HG00239", "HG00240", "HG00242", "HG00243", "HG00244", "HG00245", "HG00246", "HG00247", "HG00249", "HG00250", "HG00251", "HG00252", "HG00253", "HG00254", "HG00255", "HG00256", "HG00257", "HG00258", "HG00259", "HG00260", "HG00261", "HG00262", "HG00263", "HG00264", "HG00265", "HG00266", "HG00267", "HG00268", "HG00269", "HG00270", "HG00271", "HG00272", "HG00273", "HG00274", "HG00275", "HG00276", "HG00277", "HG00278", "HG00280", "HG00281", "HG00282", "HG00284", "HG00285", "HG00306", "HG00309", "HG00310", "HG00311", "HG00312", "HG00313", "HG00315", "HG00318", "HG00319", "HG00320", "HG00321", "HG00323", "HG00324", "HG00325", "HG00326", "HG00327", "HG00328", "HG00329", "HG00330", "HG00331", "HG00332", "HG00334", "HG00335", "HG00336", "HG00337", "HG00338", "HG00339", "HG00341", "HG00342", "HG00343", "HG00344", "HG00345", "HG00346", "HG00349", "HG00350", "HG00351", "HG00353", "HG00355", "HG00356", "HG00357", "HG00358", "HG00359", "HG00360", "HG00361", "HG00362", "HG00364", "HG00366", "HG00367", "HG00369", "HG00372", "HG00373", "HG00375", "HG00376", "HG00377", "HG00378", "HG00381", "HG00382", "HG00383", "HG00384", "HG00403", "HG00404", "HG00406", "HG00407", "HG00418", "HG00419", "HG00421", "HG00422", "HG00427", "HG00428", "HG00436", "HG00437", "HG00442", "HG00443", "HG00445", "HG00446", "HG00448", "HG00449", "HG00451", "HG00452", "HG00457", "HG00458", "HG00463", "HG00464", "HG00472", "HG00473", "HG00475", "HG00476", "HG00478", "HG00479", "HG00500", "HG00501", "HG00512", "HG00513", "HG00524", "HG00525", "HG00530", "HG00531", "HG00533", "HG00534", "HG00536", "HG00537", "HG00542", "HG00543", "HG00553", "HG00554", "HG00556", "HG00557", "HG00559", "HG00560", "HG00565", "HG00566", "HG00577", "HG00578", "HG00580", "HG00581", "HG00583", "HG00584", "HG00589", "HG00590", "HG00592", "HG00593", "HG00595", "HG00596", "HG00607", "HG00608", "HG00610", "HG00611", "HG00613", "HG00614", "HG00619", "HG00620", "HG00625", "HG00626", "HG00628", "HG00629", "HG00634", "HG00635", "HG00637", "HG00638", "HG00640", "HG00641", "HG00650", "HG00651", "HG00653", "HG00654", "HG00656", "HG00657", "HG00662", "HG00663", "HG00671", "HG00672", "HG00683", "HG00684", "HG00689", "HG00690", "HG00692", "HG00693", "HG00698", "HG00699", "HG00701", "HG00702", "HG00704", "HG00705", "HG00707", "HG00708", "HG00731", "HG00732", "HG00734", "HG00736", "HG00737", "HG00740", "HG01047", "HG01048", "HG01051", "HG01052", "HG01055", "HG01060", "HG01061", "HG01066", "HG01067", "HG01069", "HG01070", "HG01072", "HG01073", "HG01075", "HG01079", "HG01080", "HG01082", "HG01083", "HG01085", "HG01095", "HG01097", "HG01098", "HG01101", "HG01102", "HG01104", "HG01105", "HG01107", "HG01108", "HG01112", "HG01113", "HG01124", "HG01125", "HG01133", "HG01134", "HG01136", "HG01137", "HG01140", "HG01148", "HG01149", "HG01167", "HG01168", "HG01170", "HG01171", "HG01173", "HG01174", "HG01176", "HG01183", "HG01187", "HG01188", "HG01190", "HG01191", "HG01197", "HG01198", "HG01204", "HG01250", "HG01251", "HG01257", "HG01259", "HG01271", "HG01272", "HG01274", "HG01275", "HG01277", "HG01278", "HG01334", "HG01342", "HG01344", "HG01345", "HG01350", "HG01351", "HG01353", "HG01354", "HG01356", "HG01357", "HG01359", "HG01360", "HG01365", "HG01366", "HG01374", "HG01375", "HG01377", "HG01378", "HG01383", "HG01384", "HG01389", "HG01390", "HG01437", "HG01440", "HG01441", "HG01455", "HG01456", "HG01461", "HG01462", "HG01465", "HG01488", "HG01489", "HG01491", "HG01492", "HG01494", "HG01495", "HG01497", "HG01498", "HG01515", "HG01516", "HG01518", "HG01519", "HG01521", "HG01522", "HG01550", "HG01551", "HG01617", "HG01618", "HG01619", "HG01620", "HG01623", "HG01624", "HG01625", "HG01626", "NA06984", "NA06986", "NA06989", "NA06994", "NA07000", "NA07037", "NA07048", "NA07051", "NA07056", "NA07347", "NA07357", "NA10847", "NA10851", "NA11829", "NA11830", "NA11831", "NA11843", "NA11892", "NA11893", "NA11894", "NA11919", "NA11920", "NA11930", "NA11931", "NA11932", "NA11933", "NA11992", "NA11993", "NA11994", "NA11995", "NA12003", "NA12004", "NA12006", "NA12043", "NA12044", "NA12045", "NA12046", "NA12058", "NA12144", "NA12154", "NA12155", "NA12249", "NA12272", "NA12273", "NA12275", "NA12282", "NA12283", "NA12286", "NA12287", "NA12340", "NA12341", "NA12342", "NA12347", "NA12348", "NA12383", "NA12399", "NA12400", "NA12413", "NA12489", "NA12546", "NA12716", "NA12717", "NA12718", "NA12748", "NA12749", "NA12750", "NA12751", "NA12761", "NA12763", "NA12775", "NA12777", "NA12778", "NA12812", "NA12814", "NA12815", "NA12827", "NA12829", "NA12830", "NA12842", "NA12843", "NA12872", "NA12873", "NA12874", "NA12889", "NA12890", "NA18486", "NA18487", "NA18489", "NA18498", "NA18499", "NA18501", "NA18502", "NA18504", "NA18505", "NA18507", "NA18508", "NA18510", "NA18511", "NA18516", "NA18517", "NA18519", "NA18520", "NA18522", "NA18523", "NA18525", "NA18526", "NA18527", "NA18528", "NA18530", "NA18532", "NA18534", "NA18535", "NA18536", "NA18537", "NA18538", "NA18539", "NA18541", "NA18542", "NA18543", "NA18544", "NA18545", "NA18546", "NA18547", "NA18548", "NA18549", "NA18550", "NA18552", "NA18553", "NA18555", "NA18557", "NA18558", "NA18559", "NA18560", "NA18561", "NA18562", "NA18563", "NA18564", "NA18565", "NA18566", "NA18567", "NA18570", "NA18571", "NA18572", "NA18573", "NA18574", "NA18576", "NA18577", "NA18579", "NA18582", "NA18592", "NA18593", "NA18595", "NA18596", "NA18597", "NA18599", "NA18602", "NA18603", "NA18605", "NA18606", "NA18608", "NA18609", "NA18610", "NA18611", "NA18612", "NA18613", "NA18614", "NA18615", "NA18616", "NA18617", "NA18618", "NA18619", "NA18620", "NA18621", "NA18622", "NA18623", "NA18624", "NA18626", "NA18627", "NA18628", "NA18630", "NA18631", "NA18632", "NA18633", "NA18634", "NA18635", "NA18636", "NA18637", "NA18638", "NA18639", "NA18640", "NA18641", "NA18642", "NA18643", "NA18645", "NA18647", "NA18740", "NA18745", "NA18747", "NA18748", "NA18749", "NA18757", "NA18853", "NA18856", "NA18858", "NA18861", "NA18867", "NA18868", "NA18870", "NA18871", "NA18873", "NA18874", "NA18907", "NA18908", "NA18909", "NA18910", "NA18912", "NA18916", "NA18917", "NA18923", "NA18924", "NA18933", "NA18934", "NA18939", "NA18940", "NA18941", "NA18942", "NA18943", "NA18944", "NA18945", "NA18946", "NA18947", "NA18948", "NA18949", "NA18950", "NA18951", "NA18952", "NA18953", "NA18954", "NA18956", "NA18957", "NA18959", "NA18960", "NA18961", "NA18962", "NA18963", "NA18964", "NA18965", "NA18966", "NA18968", "NA18971", "NA18973", "NA18974", "NA18975", "NA18976", "NA18977", "NA18978", "NA18980", "NA18981", "NA18982", "NA18983", "NA18984", "NA18985", "NA18986", "NA18987", "NA18988", "NA18989", "NA18990", "NA18992", "NA18994", "NA18995", "NA18998", "NA18999", "NA19000", "NA19002", "NA19003", "NA19004", "NA19005", "NA19007", "NA19009", "NA19010", "NA19012", "NA19020", "NA19028", "NA19035", "NA19036", "NA19038", "NA19041", "NA19044", "NA19046", "NA19054", "NA19055", "NA19056", "NA19057", "NA19058", "NA19059", "NA19060", "NA19062", "NA19063", "NA19064", "NA19065", "NA19066", "NA19067", "NA19068", "NA19070", "NA19072", "NA19074", "NA19075", "NA19076", "NA19077", "NA19078", "NA19079", "NA19080", "NA19081", "NA19082", "NA19083", "NA19084", "NA19085", "NA19087", "NA19088", "NA19093", "NA19095", "NA19096", "NA19098", "NA19099", "NA19102", "NA19107", "NA19108", "NA19113", "NA19114", "NA19116", "NA19117", "NA19118", "NA19119", "NA19121", "NA19129", "NA19130", "NA19131", "NA19137", "NA19138", "NA19146", "NA19147", "NA19149", "NA19150", "NA19152", "NA19160", "NA19171", "NA19172", "NA19175", "NA19185", "NA19189", "NA19190", "NA19197", "NA19198", "NA19200", "NA19204", "NA19207", "NA19209", "NA19213", "NA19222", "NA19223", "NA19225", "NA19235", "NA19236", "NA19247", "NA19248", "NA19256", "NA19257", "NA19307", "NA19308", "NA19309", "NA19310", "NA19311", "NA19312", "NA19313", "NA19315", "NA19316", "NA19317", "NA19318", "NA19319", "NA19321", "NA19324", "NA19327", "NA19328", "NA19331", "NA19332", "NA19334", "NA19338", "NA19346", "NA19347", "NA19350", "NA19351", "NA19352", "NA19355", "NA19359", "NA19360", "NA19371", "NA19372", "NA19373", "NA19374", "NA19375", "NA19376", "NA19377", "NA19379", "NA19380", "NA19381", "NA19382", "NA19383", "NA19384", "NA19385", "NA19390", "NA19391", "NA19393", "NA19394", "NA19395", "NA19396", "NA19397", "NA19398", "NA19399", "NA19401", "NA19403", "NA19404", "NA19428", "NA19429", "NA19430", "NA19431", "NA19434", "NA19435", "NA19436", "NA19437", "NA19438", "NA19439", "NA19440", "NA19443", "NA19444", "NA19445", "NA19446", "NA19448", "NA19449", "NA19451", "NA19452", "NA19453", "NA19455", "NA19456", "NA19457", "NA19461", "NA19462", "NA19463", "NA19466", "NA19467", "NA19468", "NA19469", "NA19470", "NA19471", "NA19472", "NA19473", "NA19474", "NA19625", "NA19648", "NA19651", "NA19652", "NA19654", "NA19655", "NA19657", "NA19660", "NA19661", "NA19663", "NA19664", "NA19672", "NA19675", "NA19676", "NA19678", "NA19679", "NA19681", "NA19682", "NA19684", "NA19685", "NA19700", "NA19701", "NA19703", "NA19704", "NA19707", "NA19711", "NA19712", "NA19713", "NA19716", "NA19717", "NA19719", "NA19720", "NA19722", "NA19723", "NA19725", "NA19726", "NA19728", "NA19729", "NA19731", "NA19732", "NA19734", "NA19735", "NA19737", "NA19738", "NA19740", "NA19741", "NA19746", "NA19747", "NA19749", "NA19750", "NA19752", "NA19753", "NA19755", "NA19756", "NA19758", "NA19759", "NA19761", "NA19762", "NA19764", "NA19770", "NA19771", "NA19773", "NA19774", "NA19776", "NA19777", "NA19779", "NA19780", "NA19782", "NA19783", "NA19785", "NA19786", "NA19788", "NA19789", "NA19794", "NA19795", "NA19818", "NA19819", "NA19834", "NA19835", "NA19900", "NA19901", "NA19904", "NA19908", "NA19909", "NA19914", "NA19916", "NA19917", "NA19920", "NA19921", "NA19922", "NA19923", "NA19982", "NA19984", "NA19985", "NA20126", "NA20127", "NA20276", "NA20278", "NA20281", "NA20282", "NA20287", "NA20289", "NA20291", "NA20294", "NA20296", "NA20298", "NA20299", "NA20314", "NA20317", "NA20322", "NA20332", "NA20334", "NA20336", "NA20339", "NA20340", "NA20341", "NA20342", "NA20344", "NA20346", "NA20348", "NA20351", "NA20356", "NA20357", "NA20359", "NA20363", "NA20412", "NA20414", "NA20502", "NA20503", "NA20504", "NA20505", "NA20506", "NA20507", "NA20508", "NA20509", "NA20510", "NA20512", "NA20513", "NA20515", "NA20516", "NA20517", "NA20518", "NA20519", "NA20520", "NA20521", "NA20522", "NA20524", "NA20525", "NA20527", "NA20528", "NA20529", "NA20530", "NA20531", "NA20532", "NA20533", "NA20534", "NA20535", "NA20536", "NA20537", "NA20538", "NA20539", "NA20540", "NA20541", "NA20542", "NA20543", "NA20544", "NA20581", "NA20582", "NA20585", "NA20586", "NA20588", "NA20589", "NA20752", "NA20753", "NA20754", "NA20755", "NA20756", "NA20757", "NA20758", "NA20759", "NA20760", "NA20761", "NA20765", "NA20766", "NA20768", "NA20769", "NA20770", "NA20771", "NA20772", "NA20773", "NA20774", "NA20775", "NA20778", "NA20783", "NA20785", "NA20786", "NA20787", "NA20790", "NA20792", "NA20795", "NA20796", "NA20797", "NA20798", "NA20799", "NA20800", "NA20801", "NA20802", "NA20803", "NA20804", "NA20805", "NA20806", "NA20807", "NA20808", "NA20809", "NA20810", "NA20811", "NA20812", "NA20813", "NA20814", "NA20815", "NA20816", "NA20818", "NA20819", "NA20826", "NA20828"};

int main(int argc, const char * argv[])
{
    
    char  data[20000];
    
    bool male[NO_ALL];
    
    
    if(!argv[1] || !argv[2])
    {
        cout << "Input all required parameters:" << endl;
        cout << "[input-Y-VCF] [output-Y-VCF]" << endl;
        cout << "where:\n[input-Y-VCF] name of the original file with the VCF file of the Y chromosome," << endl;
        cout << "[output-Y-VCF] name of the output file with modified VCF" << endl;
        exit(1);

    }
    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
        printf("Cannot open %s\n", argv[1]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }

    FILE *out = fopen(argv[2], "w");
    
    if (out == NULL) {
        printf("Cannot open %s\n", argv[2]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }

        

    
    
    //go through the initial header
    do {
        fscanf(file, "%s", data);
    }
    while (strncmp(data, "#CHROM", 6) != 0);
    //go through all data up to names of idividuals
    fscanf(file, "%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s");
    
    size_t IndStart = my_ftell(file);
    fseek(file, 0, SEEK_SET);
    
    char * buffor = new char[IndStart];
    fread(buffor, 1, IndStart, file);
    fwrite(buffor, 1, IndStart, out);
    
    fscanf(file, "%s", data);
    for (int i = 0; i < NO_ALL; i++)
    {
        fprintf(out, "\t%s", all_ind[i]);
        if(strcmp(data, all_ind[i]) == 0)
        {
            male[i] = true;
            fscanf(file, "%s", data);
            
        }
        else
            male[i] = false;
    }
    
    fprintf(out, "\n");
    
    int formatFields = 0;
    char * ptr;

    fscanf(file, "%s", data);
    while(!feof(file))
    {
        //go through all fields up to FORMAT
        fprintf(out, "%s", data);
        for(int i = 0; i < 7; i++)
        {
            fscanf(file, "%s", data);
            fprintf(out, "\t%s", data);
        }
        //copy FORMAT field and count ':'
        fscanf(file, "%s", data);
        fprintf(out, "\t%s", data);
        ptr = strtok(data, ":");
        formatFields = 0;
        while (ptr) {
            ptr = strtok(NULL, ":");
            formatFields++;
        }
        
        
        
        //insert genotypes for all individuals
        fscanf(file, "%s", data);
        for (int i = 0; i < NO_ALL; i++)
        {
            if(male[i])
            {
                fprintf(out, "\t%s", data);
                fscanf(file, "%s", data);
            }
            else
            {
                if(formatFields == 5)
                    fprintf(out, "\t0:0:0:0:0");
                else
                {
                    fprintf(out, "\t0");
                    for(int j = 1; j < formatFields; j++)
                        fprintf(out, ":0");

                }
            }
        }
        fprintf(out, "\n");
        fscanf(file, "%s", data); //discard last (for NA21313)
        
    }
    
    
    return 0;
}

