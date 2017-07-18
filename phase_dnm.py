''' 
@queenjobo
17/07/2017

Wrapper script to use DNG phaser from trio VCF

python phase_dnm.py --vcf VCF [--probandid PROBANDID] [--motherid MOTHERID] [--fatherid FATHERID] [--parseoutput]

python /nfs/users/nfs_j/jk18/PhD/DeNovoPhaser/phase_dnm.py --vcf /lustre/scratch115/projects/ddd/users/jk18/mutatpheno/triovcfs/DDDP101263_wgs_trio.recode.vcf --dnmfile /lustre/scratch115/projects/ddd/users/jk18/mutatpheno/DNG/DDD_MAIN6529016_dnm.tab --probandbam /lustre/scratch115/projects/ddd/users/jk18/mutatpheno/bams/ROI__windel_DDDP101263_proband.bam --probandid DDD_MAIN6529016 --motherid DDD_MAIN6529046 --fatherid DDD_MAIN6529041
'''

#---------------IMPORTS---------------------

import vcf
import argparse

DNG_PATH = "/software/ddd/internal/dng-pipeline/denovogear-0.5.4_wtsi_sc13_no_limits/bin/denovogear"
#---------------FUNCTIONS-------------------

def get_options():
    '''parse command line options'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf',type = str,help = 'path to multisample VCF/ trio vcf',required = True)
    parser.add_argument('--probandbam',type = str, help = 'path to bam for proband',required = True)
    parser.add_argument('--dnmfile', type = str, help = 'dnms file', required = True)
    parser.add_argument('--probandid',type = str,help = 'sample id for proband in VCF',required = False,default = "")
    parser.add_argument('--fatherid',type = str,help = 'sample id for father in VCF',required = False,default = "")
    parser.add_argument('--motherid',type = str,help = 'sample id for mother in VCF',required = False,default = "")
    parser.add_argument('--parseoutput', help = "parse DNG output to something more user friendly", required = False,default = False, action = 'store_true')
    args = parser.parse_args()
    return(args)

def get_ids(pid,fid,mid):
    '''get ids from parsing input otherwise assume vcf is in correct order'''
    if len(pid)+len(fid)+len(mid) == 0:
        print("Assuming samples are in order proband,father,mother in VCF")
        vcf_reader = vcf.Reader(filename=myvcf,compressed=True)
        mysamples = vcf_reader.samples
        pid = mysamples[0]
        fid = mysamples[1]
        mid = mysamples[2]
    else:
        print("Leaving sample ids as entered")
    return(pid,fid,mid)

def get_geno(record,id):
    '''get genotype in form of AA'''
    call = record.genotype(id)
    bases = call.gt_bases
    if bases is not None and len(bases) == 3:
        geno = "".join(bases.split("/"))
    else:
        geno = None
    return(geno)

def make_gts(myvcf,pid,fid,mid):
    '''make gts file given vcf and ids of samples in vcf'''
    outputfile = pid + "_gts.tab"
    print("Making gts file, writing to "+outputfile)
    with open(outputfile,'w') as f:
        vcf_reader = vcf.Reader(filename=myvcf,compressed=False)
        for record in vcf_reader:
            chrm = record.CHROM
            pos = record.POS
            pgeno = get_geno(record,pid)
            fgeno = get_geno(record,fid)
            mgeno = get_geno(record,mid)
            if pgeno is None or fgeno is None or mgeno is None:
                continue
            myline = "\t".join([chrm,str(pos),pgeno,fgeno,mgeno])+"\n"
            f.write(myline)
    return(outputfile)

def run_dng_phase(probandid,dnmfile,gtsfile,probandbam):
    print("submitting dng phase")
    outfile = probandid + "_dnms_phased.dnm"
    myc = DNG_PATH + " phaser --dnm " + dnmfile + " --pgt " + gtsfile + " --bam " + probandbam + " --window 1000 > " + outfile 
    print(myc)

def main():
    args = get_options()
    pid,fid,mid = get_ids(args.probandid,args.fatherid,args.motherid)
    gtsfile = make_gts(args.vcf,pid,fid,mid)
    run_dng_phase(args.probandid,args.dnmfile,gtsfile,args.probandbam)
    
    
    

#---------------SCRIPT----------------------

if __name__ == "__main__":
    main()