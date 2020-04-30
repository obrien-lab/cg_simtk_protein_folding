#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time
import numpy as np
import mdtraj as md
import parmed as pmd

usage='Usage: python calc_cont_synth_qbb.py <domain_def> <secondary_structure_def> <eject_offset> <trun_qbb> <Q_bulk_list> <num_proc> <work_folder_list> <if_normalize>'

domain_def = None
secondary_structure_def = None
eject_offset = 30
np = 20
if len(sys.argv) == 9:
    domain_def = sys.argv[1]
    secondary_structure_def = sys.argv[2]
    eject_offset = int(sys.argv[3])
    trun_qbb_file = sys.argv[4]
    Q_bulk_list = [float(q) for q in sys.argv[5].strip().split()]
    np = int(sys.argv[6])
    work_folder_list = sys.argv[7].strip().split()
    if_normalize = int(sys.argv[8])
elif len(sys.argv) == 8:
    domain_def = sys.argv[1]
    secondary_structure_def = sys.argv[2]
    eject_offset = int(sys.argv[3])
    trun_qbb_file = sys.argv[4]
    Q_bulk_list = [float(q) for q in sys.argv[5].strip().split()]
    np = int(sys.argv[6])
    work_folder_list = ['./']
    if_normalize = int(sys.argv[7])
else:
    print("Wrong number of argv")
    print(usage)
    sys.exit()

###### parse info.log ######
def parse_info():
    num_traj = None
    prot_prefix = None
    total_chain_length = None
    fo = open('info.log')
    for line in fo:
        line = line.strip()
        if line.startswith('Number of trajectories: '):
            words = line.split()
            num_traj = int(words[-1])
        elif line.startswith('Protein psf file: '):
            words = line.split()
            prot_prefix = words[-1].split('.psf')[0]
        elif line.startswith('Total nascent chain length: '):
            words = line.split()
            total_chain_length = int(words[-1])
        elif line.startswith('File save steps: '):
            words = line.split()
            nsave = int(words[-1])
        elif line.startswith('Time step: '):
            words = line.split()
            dt = float(words[-2])
        elif line.startswith('Scaling factor: '):
            words = line.split()
            alpha = float(words[-1])
    return (num_traj, prot_prefix, total_chain_length, nsave, dt, alpha)
###### END parse info.log ######

###### parse *.out ######
def parse_out(traj_dir, nsave, dt, alpha):
    TS = []
    data_file = './output/'+traj_dir+'.out'
    traj_prefix = './traj/'+traj_dir+'/'
    init_t = 0
    tag_ter = 0
    fo = open(data_file, 'r')
    for line in fo:
        line = line.strip()
        if line.startswith('--> Elongation at length'):
            idx = int(line.split()[4])-1
            TS.append([[],[],[]])
        elif line.startswith('Simulation steps for in silico dwell time before peptidyl transfer: '):
            step = int(line.split(': ')[1])
            if step < 0:
                step = 1
            if step < nsave:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_1_final.cor']
                ts = [step*dt/1000*alpha/1e9+init_t] # in s
            elif step%nsave == 0:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_1.dcd']
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
            else:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_1.dcd', traj_prefix+'rnc_l'+str(idx+1)+'_stage_1_final.cor']
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
                ts.append(step*dt/1000*alpha/1e9+init_t)
            if idx > 0:
                init_t += step * dt / 1000 * alpha / 1e9
                TS[idx-1][2] = [traj_name, ts, [0 for i in range(len(ts))]]
        elif line.startswith('Simulation steps for in silico dwell time before translocation: '):
            step = int(line.split(': ')[1])
            if step < 0:
                step = 1
            if step < nsave:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_2_final.cor']
                ts = [step*dt/1000*alpha/1e9+init_t] # in s
            elif step%nsave == 0:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_2.dcd']
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
            else:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_2.dcd', traj_prefix+'rnc_l'+str(idx+1)+'_stage_2_final.cor']
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
                ts.append(step*dt/1000*alpha/1e9+init_t)
            init_t += step * dt / 1000 * alpha / 1e9
            TS[idx][0] = [traj_name, ts, [0 for i in range(len(ts))]]
        elif line.startswith('Simulation steps for in silico dwell time before next tRNA binding: '):
            step = int(line.split(': ')[1])
            if step < 0:
                step = 1
            if step < nsave:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_3_final.cor']
                ts = [step*dt/1000*alpha/1e9+init_t] # in s
            elif step%nsave == 0:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_3.dcd']
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
            else:
                traj_name = [traj_prefix+'rnc_l'+str(idx+1)+'_stage_3.dcd', traj_prefix+'rnc_l'+str(idx+1)+'_stage_3_final.cor']
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
                ts.append(step*dt/1000*alpha/1e9+init_t)
            init_t += step * dt / 1000 * alpha / 1e9
            TS[idx][1] = [traj_name, ts, [0 for i in range(len(ts))]]
        elif line.startswith('--> Elongation termination'):
            idx += 1
            TS.append([[],[]])
            tag_ter = 1
        elif line.startswith('Done at step ') and tag_ter == 1:
            step = int(line.split()[-1])
            if step < nsave:
                traj_name = traj_prefix+'rnc_l'+str(idx)+'_ejection_final.cor'
                ts = [step*dt/1000*alpha/1e9+init_t] # in s
            else:
                traj_name = traj_prefix+'rnc_l'+str(idx)+'_ejection.dcd'
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
            init_t += step * dt / 1000 * alpha / 1e9
            tag_ter = 2
            TS[idx][0] = [traj_name, ts, [0 for i in range(len(ts))]]
        elif line.startswith('Done at step ') and tag_ter == 2:
            step = int(line.split()[-1])
            if step < nsave:
                traj_name = traj_prefix+'rnc_l'+str(idx)+'_dissociation_final.cor'
                ts = [step*dt/1000*alpha/1e9+init_t] # in s
            else:
                traj_name = traj_prefix+'rnc_l'+str(idx)+'_dissociation.dcd'
                ts = [5000*i*dt/1000*alpha/1e9+init_t for i in range(1, int(step/nsave)+1)]
            init_t += step * dt / 1000 * alpha / 1e9
            TS[idx][1] = [traj_name, ts, [0 for i in range(len(ts))]]
    return TS
###### END parse *.out ######
    
def analysis(root_dir, ana_dir, traj_dir, init_idx, prefix, max_chain_length, eject_offset, termination_list, prot_prefix):
    global nsave, dt, alpha
    native_psf_pmd = pmd.load_file(prot_prefix+'.psf')
    native_coor = pmd.load_file(prot_prefix+'.cor').positions
    native_psf_pmd.positions = native_coor
    if eject_offset == 0:
        tmp_dir = '.temp_'+str(int(traj_dir)+init_idx-1)+'_f_'
    else:
        tmp_dir = '.temp_'+str(int(traj_dir)+init_idx-1)+'_e_'
    if if_normalize == 0:
        tmp_dir += 'r'
    else:
        tmp_dir += 'n'
    if os.path.exists(tmp_dir):
        os.system('rm -rf '+tmp_dir)
    os.mkdir(tmp_dir)
    if if_normalize == 1:
        qbb_file_name = root_dir+'/'+ana_dir+'/qbb_norm_'+str(int(traj_dir)+init_idx-1)+'_t.dat'
    else:
        qbb_file_name = root_dir+'/'+ana_dir+'/qbb_raw_'+str(int(traj_dir)+init_idx-1)+'_t.dat'
    f_qbb = open(qbb_file_name, 'w')
    f_qbb.close()
    
    TS = parse_out(traj_dir, nsave, dt, alpha)

    for chain_length in range(eject_offset+4, max_chain_length+1):
        for stage in range(3):
            if stage == 2:
                psf_name = './traj/'+traj_dir+'/'+prefix+str(chain_length+1)+'.psf'
            else:
                psf_name = './traj/'+traj_dir+'/'+prefix+str(chain_length)+'.psf'
            if chain_length == max_chain_length and stage == 2:
                break
            max_eject_idx = max([chain_length - eject_offset, 1]) - 1
            select = [i for i in range(max_eject_idx+1)]
            if if_normalize == 1:
                new_psf_pmd = native_psf_pmd[0:max_eject_idx+1]
                new_psf_pmd.save(tmp_dir+'/native.cor', format='charmmcrd', overwrite=True)
                ref_cor = tmp_dir+'/native.cor'
            else:
                ref_cor = prot_prefix+'.cor'
            qbb_list = []
            for traj_name in TS[chain_length-1][stage][0]:
                if traj_name.endswith('.cor'):
                    if os.path.exists(traj_name):
                        if not os.path.getsize(traj_name):
                            qbb_list.append('NAN')
                            continue
                    else:
                        qbb_list.append('NAN')
                        continue
                    traj_psf_pmd = pmd.load_file(psf_name)
                    traj_coor = pmd.load_file(traj_name).positions
                    traj_psf_pmd.positions = traj_coor
                    new_traj_psf_pmd = traj_psf_pmd[0:max_eject_idx+1]
                    new_traj_psf_pmd.save(tmp_dir+'/traj.cor', format='charmmcrd', overwrite=True)
                    os.system('calc_native_contact_fraction_v2.pl -i '+ref_cor+' -d '+domain_def+' -s '+secondary_structure_def+' -t '+tmp_dir+'/traj.cor > /dev/null')
                else:
                    if os.path.exists(traj_name):
                        if not os.path.getsize(traj_name):
                            for ind in range(len(TS[chain_length-1][stage][1])-1):
                                qbb_list.append('NAN')
                            continue
                    else:
                        for ind in range(len(TS[chain_length-1][stage][1])-1):
                            qbb_list.append('NAN')
                        continue
                    traj = md.load(traj_name, top=psf_name)
                    new_traj = traj.atom_slice(select)
                    new_traj.save(tmp_dir+'/traj.dcd', force_overwrite=True)
                    os.system('calc_native_contact_fraction_v2.pl -i '+ref_cor+' -d '+domain_def+' -s '+secondary_structure_def+' -t '+tmp_dir+'/traj.dcd > /dev/null')
                
                fo = open(tmp_dir+'/qbb_traj.dat', 'r')
                lines = fo.readlines()
                for i in range(1, len(lines)):
                    line = lines[i]
                    words = line.strip().split()
                    array = []
                    for q in words:
                        array.append(float(q))
                    qbb_list.append(array)
                fo.close()
            TS[chain_length-1][stage][2] = qbb_list
            if chain_length == eject_offset+4 and stage == 0:
                header = lines[0].strip().split()
                f_qbb = open(qbb_file_name, 'a')
                f_qbb.write('%12s %5s %5s '%('Time (s)', 'Len', 'Stage'))
                for s in header:
                    f_qbb.write('%7s '%s)
                f_qbb.write('\n')
                for i in range(eject_offset+3):
                    for j in range(len(TS[i])):
                        for k in range(len(TS[i][j][1])):
                            f_qbb.write('%12.6f %5d %5d '%(TS[i][j][1][k], i+1, j+1))
                            for m in range(len(header)):
                                f_qbb.write('%7.4f '%0)
                            f_qbb.write('\n')
                f_qbb.close()
            
            f_qbb = open(qbb_file_name, 'a')
            for k in range(len(qbb_list)):
                if qbb_list[k] == 'NAN':
                    continue
                f_qbb.write('%12.6f %5d %5d '%(TS[chain_length-1][stage][1][k], chain_length, stage+1))
                for m in range(len(qbb_list[k])):
                    f_qbb.write('%7.4f '%qbb_list[k][m])
                f_qbb.write('\n')
            f_qbb.close()
       
    select = [i for i in range(max_chain_length)]
    psf_name = './traj/'+traj_dir+'/'+prefix+str(max_chain_length)+'.psf'
    if if_normalize == 1:
        new_psf_pmd = native_psf_pmd[0:max_chain_length]
        new_psf_pmd.save(tmp_dir+'/native.cor', format='charmmcrd', overwrite=True)
        ref_cor = tmp_dir+'/native.cor'
    else:
        ref_cor = prot_prefix+'.cor'
    for stage in range(2):
        traj_name = TS[max_chain_length][stage][0]
        if traj_name.endswith('.cor'):
            traj_psf_pmd = pmd.load_file(psf_name)
            traj_coor = pmd.load_file(traj_name).positions
            traj_psf_pmd.positions = traj_coor
            new_traj_psf_pmd = traj_psf_pmd[0:max_chain_length]
            new_traj_psf_pmd.save(tmp_dir+'/traj.cor', format='charmmcrd', overwrite=True)
            os.system('calc_native_contact_fraction_v2.pl -i '+ref_cor+' -d '+domain_def+' -s '+secondary_structure_def+' -t '+tmp_dir+'/traj.cor > /dev/null')
        else:
            traj = md.load(traj_name, top=psf_name)
            new_traj = traj.atom_slice(select)
            new_traj.save(tmp_dir+'/traj.dcd', force_overwrite=True)
            os.system('calc_native_contact_fraction_v2.pl -i '+ref_cor+' -d '+domain_def+' -s '+secondary_structure_def+' -t '+tmp_dir+'/traj.dcd > /dev/null')
        fo = open(tmp_dir+'/qbb_traj.dat', 'r')
        lines = fo.readlines()
        qbb_list = []
        for i in range(1, len(lines)):
            line = lines[i]
            words = line.strip().split()
            array = []
            for q in words:
                array.append(float(q))
            qbb_list.append(array)
        fo.close()
        TS[max_chain_length][stage][2] = qbb_list
        f_qbb = open(qbb_file_name, 'a')
        for k in range(len(qbb_list)):
            f_qbb.write('%12.6f %5d %5d '%(TS[max_chain_length][stage][1][k], max_chain_length, stage+4))
            for m in range(len(qbb_list[k])):
                f_qbb.write('%7.4f '%qbb_list[k][m])
            f_qbb.write('\n')
        f_qbb.close()
    
    os.system('rm -rf '+tmp_dir)
    
    

###################### MAIN ################################
if len(work_folder_list) == 0:
    work_folder_list = ['./']
root_dir = os.getcwd()
if eject_offset == 0:
	ana_dir = 'qbb_full_vs_T'
else:
	ana_dir = 'qbb_eject_vs_T'
if not os.path.exists(ana_dir):
    os.mkdir(ana_dir)

max_length = 0
max_num_traj = 0
for work_folder in work_folder_list:
    os.chdir(work_folder)
    if work_folder == './':
        init_idx = 1
    else:
        words = work_folder.split('-')
        init_idx = int(words[0])
    (num_traj, prot_prefix, total_chain_length, nsave, dt, alpha) = parse_info()

    prefix = 'rnc_l'

    traj_dir_list = []
    for file in os.listdir('./traj/'):
        if os.path.isdir('./traj/'+file):
            traj_dir_list.append(file)    
    
    pool = multiprocessing.Pool(np)
    for traj_dir in traj_dir_list:
        chain_length_list = []
        termination_list = []
        for file in os.listdir('./traj/'+traj_dir):
            if file.startswith(prefix) and file.endswith('_stage_3_final.cor'):
                name = os.path.splitext(file)[0]
                words = name.split('l')
                words = words[-2].split('_')
                chain_length_list.append(int(words[0]))
            elif file.startswith(prefix) and (file.endswith('_ejection.dcd') or file.endswith('_dissociation.dcd')):
                termination_list.append(file)
        max_chain_length = max(chain_length_list)
        if max_chain_length > max_length:
            max_length = max_chain_length
        max_num_traj += 1
        pool.apply_async(analysis, (root_dir, ana_dir, traj_dir, init_idx, prefix, max_chain_length, eject_offset, termination_list, prot_prefix))

    while True:
        time.sleep(10)
        if len(pool._cache) == 0:
            print('All Done.')
            break

    pool.close()
    pool.join()

    os.chdir(root_dir)

#os.chdir(root_dir+'/'+ana_dir+'/')
#fo = open('plot_qbb.m', 'w')
#fo.write('''clear all
#dt = '''+str(dt)+'''/1000;
#nsave = '''+str(nsave)+''';
#alpha = '''+str(alpha)+''';
#dt = dt*nsave*alpha/1e6;
#n_traj = '''+str(max_num_traj)+''';
#max_length = '''+str(max_length)+''';
#Q_bulk_list = [''');
#for i in Q_bulk_list:
#    fo.write(str(i)+', ')
#fo.write('''];
#Q_length_file = \''''+trun_qbb_file+'''\';
#eject_offset = '''+str(eject_offset)+''';
#domain_idx = find(Q_bulk_list ~= 1);
#D = {};
#D_ter = {};
#for i=1:n_traj
#C = importdata(['qbb_cs_' num2str(i) '.dat'], ' ', 1);
#D{i} = C.data(:,2:end);
#domain_str = C.colheaders;
#C = importdata(['qbb_te_' num2str(i) '.dat'], ' ', 1);
#D_ter{i} = C.data(:,1:end);
#end
#
#if ~isequal(Q_length_file, '')
#Q_length = importdata(Q_length_file, ' ', 1);
#Q_length = Q_length.data;
#if Q_length(1,1) ~= 1
#    L = [1];
#    for i = 2:size(Q_length,2)
#        L(1,i) = 0;
#    end
#    Q_length = [L; Q_length];
#end
#if Q_length(end,1) ~= max_length
#    L = [max_length];
#    L(1,2:2:size(Q_length,2)) = Q_bulk_list;
#    L(1,3:2:size(Q_length,2)) = Q_length(end, 3:2:size(Q_length,2));
#    Q_length = [Q_length; L];
#end
#Q_avg_length = Q_length(:, 2:2:size(Q_length,2));
#Q_std_length = Q_length(:, 3:2:size(Q_length,2));
#Q_avg_length_fit = interp1(Q_length(:,1), Q_avg_length, (1:max_length)', 'pchip');
#Q_ustd_length_fit = interp1(Q_length(:,1), Q_avg_length+Q_std_length, (1:max_length)', 'pchip');
#Q_lstd_length_fit = interp1(Q_length(:,1), Q_avg_length-Q_std_length, (1:max_length)', 'pchip');
#end
#
#sub_length = 20;
#sub_height = 4;
#f_length = sub_length;
#f_height = sub_height*(length(domain_idx));
#figure(1)
#set(gcf,'Resize', 'off', 'Units','centimeters',...
#'Position',[8 1.5 f_length f_height],'paperpositionmode','auto');
#for j = 1:length(domain_idx)
#subplot(length(domain_idx), 1, j)
#hold on
#Q = [];
#max_frame = 0;
#for i=1:n_traj
#    plot(D{i}(:,domain_idx(j)), '-', 'Color', [135 206 250]/255, 'Linewidth', 1.0)
#    plot(max_length+0.5*(1:size(D_ter{i},1)), D_ter{i}(:,domain_idx(j)), '-', 'Color', [135 206 250]/255, 'Linewidth', 1.0)
#    if max_frame < size(D_ter{i},1)
#        max_frame = size(D_ter{i},1);
#    end
#end
#for k=1:max_length
#    a=0;
#    n=0;
#    for i=1:n_traj
#        if k <= length(D{i}(:,domain_idx(j)))
#            a=a+D{i}(k,domain_idx(j));
#            n=n+1;
#        end
#    end
#    Q(k)=a/n;
#end
#plot(Q,'-b','Linewidth', 1.0)
#if ~isequal(Q_length_file, '')
#    plot(eject_offset+1:max_length, Q_avg_length_fit(1:max_length-eject_offset,domain_idx(j))', '-r', 'Linewidth', 1)
#    patch([eject_offset+1:max_length max_length:-1:eject_offset+1], [Q_ustd_length_fit(1:max_length-eject_offset,domain_idx(j))' Q_lstd_length_fit(max_length-eject_offset:-1:1,domain_idx(j))'], 'r','FaceAlpha',0.2,'EdgeColor','none','LineWidth',1.0)
#    plot([max_length max_length+0.5*max_frame], [Q_bulk_list(domain_idx(j)) Q_bulk_list(domain_idx(j))],'-r','Linewidth', 1)
#else
#    plot([1 max_length+0.5*max_frame], [Q_bulk_list(domain_idx(j)) Q_bulk_list(domain_idx(j))],'-r','Linewidth', 1)
#end
#plot([max_length max_length], [0 1],'-k','Linewidth', 1.5)
#set(gca, 'fontsize',11,'fontweight','normal','LineWidth',1.5,'fontname','Nimbus Roman No9 L')
#box on
#grid on
#xlh = xlabel(' ','fontsize',13,'color','k','Interpreter','latex');
#ylh = ylabel('$Q_{\\rm syn}$','fontsize',13,'color','k','Interpreter','latex');
#text('String', 'Chain length (aa)', 'Position', [0.5*(1+max_length),1.2*xlh.Position(2)],...
#    'fontsize',13,'color','k','Interpreter','latex',...
#    'HorizontalAlignment','center','VerticalAlignment','top')
#text('String', 'Experimental time (ms)', 'Position', [0.5*(2*max_length+0.5*max_frame),...
#    1.2*xlh.Position(2)],'fontsize',13,'color','k','Interpreter','latex',...
#    'HorizontalAlignment','center','VerticalAlignment','top')
#text('String', '$Q_{\\rm ter}$', 'Position', [-0.5*ylh.Position(1)+max_length+0.5*max_frame,0.5],...
#    'fontsize',13,'color','k','Interpreter','latex',...
#    'HorizontalAlignment','center','rotation',90)
#title_str = domain_str{domain_idx(j)+1};
#if ~isempty(strfind(title_str, 'D'))
#    title_str = ['Domain ' title_str(3:end)];
#elseif ~isempty(strfind(title_str, '|'))
#    idx = strfind(title_str, '|');
#    title_str = ['Interface ' title_str(1:idx(1)-1) '$\\mid$' title_str(idx(1)+1:end)];
#else
#    title_str = [upper(title_str(1)) title_str(2:end)];
#end
#title(title_str, 'Fontsize', 13, 'Interpreter','latex')
#axis([1 max_length+0.5*max_frame 0 1])
#xtick_1 = 0:20:max_length; 
#xtick_2 = max_length+0.5*(20/dt:20/dt:max_frame);
#xtick = [xtick_1 xtick_2];
#xticklabel = {};
#for i = 1:length(xtick_1)
#    xticklabel{i} = num2str(xtick_1(i));
#end
#for i = 1:length(xtick_2)
#    xticklabel{length(xtick_1)+i} = 20*i;
#end
#set(gca, 'XTick', [xtick_1 xtick_2], 'XTickLabel', xticklabel);
#end
#saveas(gcf, 'qbb_vs_length.svg')
#quit\n''');
#fo.close()
#os.system('matlab -nodisplay -r plot_qbb > /dev/null')
