######
# Convert merged bed file to gtf file
# We need this to call peaks from 
######
import getopt, sys, os.path


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=","output="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    infile = None
    outfile = None

    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            if os.path.isfile(a):
                infile = a
        elif o in ("-o", "--output"):
            outfile = a
        else:
            assert False, "Unhandled option"

    if infile is not None and outfile is not None:
        run(infile, outfile)
    else:
        usage()


def usage():
    print("\nUsage: python bed2gtf [options] <mandatory>")
    print("Options:")
    print("\t-h, --help:\n\t\t show this help message and exit")
    print("Mandatory:")
    print("\t-i, --input:\n\t\t File with the regions in bed format")
    print("\t-o, --output:\n\t\t Name of the gtf file output file. Directory where the file will be created should exist!")

def run(infile, outfile):

    inf  = open(infile, 'r')
    outf = open(outfile,'w')

    cont = 1
    for linea in inf:
        if linea.startswith('#'):
            continue
        # print(linea)
        linea_split = linea.split()
        chrom = str(linea_split[1]).strip()
        ini_pos = str(linea_split[2]).strip()
        fin_pos = str(linea_split[3]).strip()
        strand = str(linea_split[4]).strip()        

        outf.write(f'{chrom}\tmerged\tpeak\t{ini_pos}\t{fin_pos}\t.\t{strand}\t.\tpeak_id "{chrom}_{ini_pos}_{fin_pos}_{strand}";\n')

        cont += 1

    inf.close()
    outf.close()


if __name__ == "__main__":
    main()