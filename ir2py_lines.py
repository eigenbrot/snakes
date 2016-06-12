import numpy as np
import re

def parse_fitprofs(irfile, numlines = 3):

    with open(irfile, 'r') as f:
        lines = f.readlines()
    r = re.compile('Ap\s(\d+):')
        
    aps = []
    l = 0
    while l < len(lines) - 1:
        # print l, lines[l]
        # raw_input('wait')
        if lines[l][0] == '#':
            if 'Ap' in lines[l]:
                apnum = int(r.search(lines[l]).groups()[0])
                aps.append(apnum)
                tc = []
                tw = []
                te = []
                tf = []
                tfe = []
                tb = []
                tp = []
                l += 3
                count = 0
                
            else:
                continue
                
        while count < numlines:
            vals = lines[l].split()
            tc.append(float(vals[0]))
            tw.append(float(vals[5]))
            tf.append(float(vals[2]))
            tb.append(float(vals[1]))
            tp.append(float(vals[4]))
            try:
                if '(' in lines[l+1]:
                    vals = lines[l+1].split(')')
                    te.append(float(vals[0].split('(')[1]))
                    tfe.append(float(vals[2].split('(')[1]))
                    l += 2
                elif lines[l+1][0] == '#':
                    l += 1
                else:
                    te = [99.] * numlines
                    tfe = [99.] * numlines
                    l += 1
            except IndexError:
                l += 1
            count += 1

        try:
            cents = np.vstack((cents, tc))
            widts = np.vstack((widts, tw))
            errs = np.vstack((errs, te))
            fluxes = np.vstack((fluxes, tf))
            fluxerr = np.vstack((fluxerr, tfe))
            backgrd = np.vstack((backgrd, tb))
            peaks = np.vstack((peaks, tp))
        except NameError:
            cents = np.array(tc)
            widts = np.array(tw)                
            errs = np.array(te)
            fluxes = np.array(tf)
            fluxerr = np.array(tfe)
            backgrd = np.array(tb)
            peaks = np.array(tp)

    return aps, cents, widts, errs, fluxes, fluxerr, backgrd, peaks

def compute_widths(irfile, output, numlines = 3):

    dat = parse_fitprofs(irfile, numlines = numlines)
    
    mcent = np.mean(dat[1],axis=0)
    scent = np.std(dat[1],axis=0)
    mwidt = np.mean(dat[2],axis=0)
    swidt = np.std(dat[2],axis=0)

    with open(output,'w') as f:
        f.write('# {:>8}{:>10}{:>10}{:>10}\n#\n'.\
                format('cent', 'sigcent', 'width', 'sigwidt'))
        for i in range(mcent.size):
            f.write(str('{:10.3f}'*4+'\n').\
                    format(mcent[i],scent[i],mwidt[i],swidt[i]))

    return

def batch_lines(output):

    with open(output,'w') as f:
        f.write(str('# {:>8}'+'{:>10}'*10).format('Wave',
                                                  "6''",
                                                  "5''",
                                                  "4''",
                                                  "3''",
                                                  "2''",
                                                  "sig 6''",
                                                  "sig 5''",
                                                  "sig 4''",
                                                  "sig 3''",
                                                  "sig 2''"))
        f.write('\n#\n')
        for l in range(15):
            cent_list = []
            widt_list = []
            sig_list = []
        
            for s in range(5):
                ifile = 'fitprofs/l{}_{}.fitp'.format(l+1,s+2)
                print ifile
                if l == 0:
                    dat = parse_fitprofs(ifile,4)
                else:
                    dat = parse_fitprofs(ifile,1)
            
                cent_list.append(np.mean(dat[1],axis=0))
                widt_list.append(np.mean(dat[2],axis=0))
                sig_list.append(np.std(dat[2],axis=0))
            
            if l == 0:
                cent_stack = np.vstack(cent_list).T
                widt_stack = np.vstack(widt_list).T
                sig_stack = np.vstack(sig_list).T

                for c in range(cent_stack.shape[0]):
                    cent = np.mean(cent_stack[c])
                    cetnstd = np.std(cent_stack[c])
                    args = np.r_[cent,widt_stack[c][::-1],sig_stack[c][::-1]]
                    f.write(str('{:10.3f}'*11+'\n').format(*args))
            
            else:
                cent_stack = np.hstack(cent_list)
                widt_stack = np.hstack(widt_list)
                sig_stack = np.hstack(sig_list)

                cent = np.mean(cent_stack)
                cetnstd = np.std(cent_stack)
                args = np.r_[cent,widt_stack[::-1],sig_stack[::-1]]
                f.write(str('{:10.3f}'*11+'\n').format(*args))
    return
               
