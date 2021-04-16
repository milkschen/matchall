from pysam import VariantFile

fn_panel_vcf = "/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/reference_flow/resources/1kg_vcf/annotated/chr21.vcf.gz"
fn_happy_vcf = "/scratch/groups/blangme2/naechyun/imputation_dv/deepvariant_af/HG003.deepvariant_af.happy.vcf.gz"

f_happy_vcf = VariantFile(fn_happy_vcf)

cycle = 0
outcome_tmp = {}
outcome = {
    'SNP':{
        'TRUTH.TOTAL': 0,
        'TRUTH.TP': 0,
        'TRUTH.FN': 0,
        'QUERY.TOTAL': 0,
        'QUERY.FP': 0,
        'QUERY.UNK': 0,
        'FP.gt': 0
    },
    'INDEL':{
        'TRUTH.TOTAL': 0,
        'TRUTH.TP': 0,
        'TRUTH.FN': 0,
        'QUERY.TOTAL': 0,
        'QUERY.FP': 0,
        'QUERY.UNK': 0,
        'FP.gt': 0
    }
}

def update_happy_outcome(bd, bvt, outcome):
    # TRUTH
    if bd[0] == 'TP':
        outcome[bvt[0]]['TRUTH.TP'] += 1
    elif bd[0] == 'FN':
        outcome[bvt[0]]['TRUTH.FN'] += 1
    # QUERY
    if bd[1] == 'FP':
        outcome[bvt[1]]['QUERY.FP'] += 1


for var in f_happy_vcf.fetch():
    if var.info.get('Regions'):
        bvt = []
        bd = []
        for sample in var.samples.items():
            for v_format in sample[1].items():
                if v_format[0] == 'BD':
                    bd.append(v_format[1])
                elif v_format[0] == 'BVT':
                    bvt.append(v_format[1])
            #print(f'sample name: {sample[0]}')
            #print(f'format values: {sample[1].items()}')
        update_happy_outcome(bd, bvt, outcome)
        bd = ' '.join(bd)
        bvt = ' '.join(bvt)
        for v_type in bvt:
            if outcome_tmp.get(v_type):
                if outcome_tmp[bvt].get(bd):
                    outcome_tmp[bvt][bd] += 1
                else:
                    outcome_tmp[bvt][bd] = 1
            else:
                outcome_tmp[bvt] = {bd: 1}
            # stop update if bvt[0] == bvt[1]
            if all(t == v_type for t in bvt):
                break
        #print(bd)
        if cycle > 0 and cycle % 1000000 == 0:
            #cycle = 0
            #input()
            print(f'Process {cycle} records')
        cycle += 1
print(outcome)
print(outcome_tmp)
