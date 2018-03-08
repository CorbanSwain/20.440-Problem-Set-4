import pymol

prot_seq = 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVN'
#           123456789012345678901234567890123456789012345678901234567
#           0        1         2         3         4         5

views = [(-0.361983389,    0.388416916,   -0.847408175,
           0.437206864,    0.873611212,    0.213667795,
           0.823295772,   -0.293148011,   -0.486049443,
           0.000000000,    0.000000000,  -85.331138611,
           3.321181297,   26.772531509,    9.855480194,
           28.665481567,  141.996810913,  -20.000000000)]

debug_mode = False;
blocks = [(4, 11), (16, 28), (35, 39), (44, 53)]
fully_conserved = [17, 22, 36, 46, 49]
prot_sel = '(c. A)'
ligand_sel = '(c. B)'

def setup():
    cmd.reset() # reset view
    hide_all()
    util.performance(0) # maximum quality
    cmd.bg_color('white') # white background
    cmd.set('cartoon_fancy_helices','1',quiet=0) # fancy helixes
    cmd.set('surface_color_smoothing_threshold', 0.7)
    cmd.set('surface_color_smoothing',20)
    cmd.set('transparency_mode',3,quiet=0)
    cmd.set('surface_quality',2)
    cmd.set('surface_cavity_mode', 2)
    cmd.set('cavity_cull', 150)

def hide_all():
    cmd.hide('everything', '*')
    
def color_blocks():
    hide_all()
    main_color = 'cyan'
    ligand_color = 'gray90'
    color_list = ['red', 'forest', 'orange', 'magenta']
    cmd.color(main_color, prot_sel)
    cmd.color(ligand_color, ligand_sel)
    cmd.show('cartoon', prot_sel)
    cmd.show('cartoon', ligand_sel)
    for i, (iStart, iEnd) in enumerate(blocks):
        if debug_mode:
            print('%d - %s - %s' % (i, resi_sel(iStart, iEnd), color_list[i]))
        cmd.color(color_list[i], resi_sel(iStart, iEnd))

def show_consv():
    for resi in fully_conserved:
        cmd.show('stick', '%s & i. %d & !(n. C+N+O) & !(e. H)' % (prot_sel, resi))
        if prot_seq[resi - 1] == 'G':
            cmd.show('stick', '%s & i. %d & e. H' % (prot_sel, resi))
        if prot_seq[resi - 1] == 'P':
            cmd.show('stick', '%s & i. %d & n. N'  % (prot_sel, resi))

def show_view(num):
    cmd.set_view(views[num])
        
def resi_sel(start, end):
    nums = range(start, end + 1)
    nums_str = '+'.join(['%d' % n for n in nums])
    return '%s &  i. %s' % (prot_sel, nums_str)

def save_im(name):
    cmd.set('ray_trace_mode', 3)
    cmd.set('ray_texture', 1)
    cmd.ray(1800, 1200)
    cmd.png('protein_figure_' + name);
    print('Saving Done')
    
