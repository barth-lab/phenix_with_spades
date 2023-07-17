
"""
$PHENIX/phenix/phenix/rosetta/utils.py
Part of phenix autobuild. Utilities for running Rosetta from Phenix
2010-08-04
TT
"""
from __future__ import division
from __future__ import print_function

from libtbx.utils import Sorry
import libtbx.callbacks # import dependency
import sys, os
from phenix.autosol.init_utils import init_utils
from phenix.autosol.get_pdb_inp import get_pdb_inp
import iotbx.pdb
from libtbx.utils import null_out

class utils(init_utils):

  def set_up_reflections_and_r_free(self,params,out=sys.stdout,
     require_data=True,require_data_or_map_coeffs=True):
    # convert to mtz; define free set if necessary;
    print("\n",80*"=", file=out)
    print("Setting up reflection file and labels", file=out)
    print(80*"=", file=out)

    if hasattr(params,'refinement') and \
       getattr(params.refinement,'use_hl',None) == False and \
       hasattr(params.input_files,'hl_labels') and \
       params.input_files.hl_labels:
      raise Sorry("Sorry you cannot specify hl_labels if use_hl=False")

    if hasattr(params.crystal_info,'correct_aniso') and \
         params.crystal_info.correct_aniso:
      print("\nApplying anisotropy correction to data", file=out)
      params.crystal_info.correct_aniso=False
      self.set_up_reflections_and_r_free(params,out=out,
         require_data=require_data,
         require_data_or_map_coeffs=require_data_or_map_coeffs)
      # now we have read in data in usual way. Apply correction to the
      # data and then read in again
      temp_data=os.path.join(params.directories.temp_dir,"aniso_temp_file.mtz")
      args=[params.input_files.data]
      args.append('hklout=%s' %(temp_data))
      args.append('hklout_type=mtz')
      args.append('aniso.action=remove_aniso')
      if params.crystal_info.final_b is None:
         args.append("final_b=eigen_min")
      else:
         args.append("final_b=user_b_iso")
         args.append("b_iso=%s" %(str(params.crystal_info.final_b)))
      from mmtbx.command_line import massage_data
      massage_data.run(args=args, out=null_out())
      self.wait_for_file_to_appear(temp_data,
        max_ticks=params.control.max_wait_time,
        allow_failure=False)
      # now use this for data and the original file for free R:
      params.input_files.free_r_data=params.input_files.data
      params.input_files.free_r_labels=params.input_files.data_labels # only freer used
      params.input_files.data=temp_data
      params.input_files.data_labels=None
      print("\nTemporary aniso-corrected data in %s"%(temp_data), file=out)


    from phenix.autosol.read_hkl import read_hkl
    if not params.input_files.data or \
         not os.path.isfile(params.input_files.data):
      if require_data:
        if params.input_files.data:
          raise Sorry("The data file '%s' seems to be missing " %(params.input_files.data))
        else:
          raise Sorry("Please specify a data file with a command such as "+
            "'data=data.mtz'")

      elif require_data_or_map_coeffs:
         if (not params.input_files.map_coeffs or \
              not os.path.isfile(params.input_files.map_coeffs)):
           raise Sorry("Sorry need a data file or map_coeffs file")

    # if...
    #  (1) data has FP SIGFP FreeR_flag, optionally labin is given to define
    #       no data_labels or FreeR_labels allowed, no hl_labels:
    #       -> use this as is
    #  otherwise...
    #  (2) data has F,SIGF or I,SIGI, data or  free_r_data has FreeR_flags
    #       optionally data_labels and/or free_r_labels and/or hl_labels,
    #       no labin allowed
    #       -> convert this to mtz F SIGF FreeR_flags, specify labin and data,
    #          remove data_labels and free_r_labels
    #  otherwise...
    #  (3) stop with error (say what is missing: data, freeR_flags)

    # So: if data_labels or free_r_labels set or free_r_data set
    # or hl_labels or hl_data set
    # or we cannot interpret data file then import the data:
    # Also if data has HL data and use_hl=False, import it to forget hl coeffs

    r=None
    need_to_import=False
    if params.input_files.free_r_data and \
       params.input_files.free_r_data != params.input_files.data:
      print("Free R data present and not same as data", file=out)
      need_to_import=True

    elif hasattr(params.input_files,'hl_data') and \
       params.input_files.hl_data and \
       params.input_files.hl_data != params.input_files.data:
      print("HL data present and not same as data", file=out)
      need_to_import=True

    if not need_to_import:
      r=read_hkl(params.input_files.data,split=True,quiet=True)
      if r.get_file_type() != 'ccp4_mtz':
        need_to_import=True
        print("File type is not mtz", file=out)

    if r and not need_to_import:
      self.get_free_r_and_fp_labels(r=r,out=out,require_phib=False,
          check_use_hl_anom=True,nothing_required=True)
      if self.data_fp_label in ['None',None] or \
          self.data_sigfp_label in ['None',None] or \
          self.data_freer_label in ['None',None] or \
          self.phenix_refine_label in ['None',None]:
        print("Missing labels", file=out)
        need_to_import=True
      elif params.input_files.data_labels and \
           params.input_files.data_labels != self.phenix_refine_label:
        print("Importing as target phenix_refine_label is '%s' " %(
        self.phenix_refine_label) +\
        " and data_labels are '%s'" %(params.input_files.data_labels), file=out)
        need_to_import=True
      elif params.input_files.free_r_labels and \
           params.input_files.free_r_labels != self.data_freer_label:
        print("Free R labels do not match requested", file=out)
        need_to_import=True
      elif hasattr(params.input_files,'hl_labels') and \
           params.input_files.hl_labels and \
           params.input_files.hl_labels != self.hl_label:
        print("HL labels do not match requested", file=out)
        need_to_import=True
      elif self.hl_label and hasattr(params,'refinement') \
        and getattr(params.refinement,'use_hl',None) == False:
        print("Skipping HL information", file=out)
        need_to_import=True

    if need_to_import:
      mtz_out=os.path.join(params.directories.temp_dir,
         'imported_fp_sigfp_freer.mtz')
      use_hl=None
      if hasattr(params,'refinement'):
        use_hl=getattr(params.refinement,'use_hl',None)
      if hasattr(params.input_files,'hl_data'):
        hl_data=params.input_files.hl_data
        hl_labels=params.input_files.hl_labels
      else:
        hl_data=None
        hl_labels=None
      params.input_files.data=self.import_data_and_free_r(
        data=params.input_files.data,
        data_labels=params.input_files.data_labels,
        free_r_data=params.input_files.free_r_data,
        free_r_labels=params.input_files.free_r_labels,
        hl_data=hl_data,
        hl_labels=hl_labels,
        use_hl=use_hl,
        mtz_out=mtz_out,out=out)
      params.input_files.free_r_data=None
      params.input_files.data_labels=None
      params.input_files.free_r_labels=None
      if hasattr(params.input_files,'hl_data'):
        params.input_files.hl_labels=None
      r=read_hkl(params.input_files.data,split=True,quiet=True)

    # Now we should be ready:
    self.space_group_from_data_file=r.sg # 2011-11-02 (just the symbol)
    self.crystal_symmetry=r.crystal_symmetry
    self.check_use_hl_anom(input_hkl=r,check=True)
    self.get_free_r_and_fp_labels(r=r,out=out,require_phib=False,
      check_use_hl_anom=True)
    print("LABIN LINE TO BE USED: FP=%s SIGFP=%s FreeR_flag=%s " %(
      self.data_fp_label, self.data_sigfp_label, self.data_freer_label), end=' ', file=out)
    # check to see if we are going to use hl coeffs (and what kind, if so):
    # self.hl_label or self.hl_anom_labels_from_list
    self.hl_label_to_use=None
    if hasattr(params,'refinement') and \
       getattr(params.refinement,'use_hl',None) != False:
      if getattr(self,'hl_anom_list',None) and \
          getattr(params.refinement,'use_hl_anom_in_refinement',None):
        hl_labels=self.hl_anom_list
      elif self.hl_label:
        hl_labels=self.hl_label.split(",")
      else:
        hl_labels=None
      if hl_labels:
        print("HLA=%s HLB=%s HLC=%s HLD=%s " % tuple(hl_labels), file=out)
        self.hl_label_to_use=",".join(hl_labels)
      else:
        print(file=out)
    else:
     print(file=out)

    if not self.check_for_required():
        raise Sorry("\nSorry, the data label %s" %(self.missing_label) +
         " needs to be present in the "+
            " \ninput file %s " %(params.input_files.data))


    # now read map_coeffs file if present:
    if params.input_files.map_coeffs:
      r=read_hkl(params.input_files.map_coeffs,split=True,quiet=True)
      if r.get_file_type() != 'ccp4_mtz':
        raise Sorry("\nSorry, the map coeffs file '%s' must be an mtz file\n"
           %(params.input_files.map_coeffs) )
      # get the labels
      self.get_map_labels(r,out=out)

    if hasattr(params,'crystal_info') and \
         hasattr(params.crystal_info,'ncs_copies'):
      sg=r.sg
      cell=r.cell_params
      # get possible values for ncs_copies, if not already defined
      if self.is_auto(params.crystal_info.ncs_copies) or \
          not params.crystal_info.ncs_copies:
        print("\nIdentifying plausible values of ncs_copies\n", file=out)
        possible_ncs_copies=self.guess_values_of_ncs_copies(params,
          sg,cell,out=out)
        params.crystal_info.ncs_copies=possible_ncs_copies

  def check_for_required(self):
    for item,name in zip([self.data_fp_label, self.data_sigfp_label,
         self.data_freer_label],['FP','SIGFP','FreeR_flag']):
      if str(item) == 'None':
        self.missing_label=name
        return False
    return True
  def get_fp_mtz_labels(self,labin_line):
    # purpose: return "myFP,mySIGFP" if given  "FP=myFP SIGFP=mySIGFP"
    fp_value=None
    sigfp_value=None
    for item in labin_line.split():
      spl=item.split("=")
      if len(spl)<2: continue
      if spl[0].upper()=="FP":
        fp_value=spl[1]
      elif spl[0].upper()=="SIGFP":
        sigfp_value=spl[1]
    if not fp_value or not sigfp_value:
     return None
    return " ".join([fp_value,sigfp_value])

  def refine_mr_models(self,params,mr_model_list,out=sys.stdout,
      ncs_object_list=None,
      refine_in_temp_dir=False):
    refined_mr_model_list=[]
    refined_mr_model_mtz_list=[]
    refined_labin_list=[]
    data_list=[]
    crystal_symmetry_list=[]
    labin_map_coeffs="FP=2FOFCWT PHIB=PH2FOFCWT"
    if params.non_user_params.dummy_refinement or not \
         params.place_model.refine_after_mr:
      dummy_refinement=True
    else:
      dummy_refinement=False

    if ncs_object_list is None:
       ncs_object_list=len(mr_model_list)*[None]
    for mr_model,ncs_object in zip(mr_model_list,ncs_object_list):
      pdb_input = get_pdb_inp(file_name=mr_model)
      hierarchy = pdb_input.construct_hierarchy()
      crystal_symmetry=self.none_if_no_space_group_info(
         pdb_input.crystal_symmetry_from_cryst1())
      if crystal_symmetry is None:
         raise Sorry("\nPlease add a CRYST1 card with symmetry to %s" %(
            mr_model))
      crystal_symmetry_list.append(crystal_symmetry)
      space_group_symbol=self.get_space_group_symbol(
         crystal_symmetry=crystal_symmetry)

      # refine this model
      if refine_in_temp_dir:
        from phenix.autosol.trim_file_name import trim_file_name
        new_data=os.path.join(params.directories.temp_dir,
          trim_file_name(mr_model[:-4]+"_data.mtz").trimmed_file)
        refined_prefix=os.path.join(params.directories.temp_dir,
          trim_file_name(mr_model[:-4]+"_ref").trimmed_file)
      else:
        new_data=mr_model[:-4]+"_data.mtz"
        refined_prefix=mr_model[:-4]+"_ref"
      if ncs_object is None:
        ncs=False
      else:
        ncs=True

      if self.params.control.real_space_optimize: # don't do anything
        refined_mr_model,map_coeffs=mr_model,None
      else:  # usual
        refined_mr_model,map_coeffs=self.run_refine(
         mr_model,refined_prefix,out=out,
         crystal_symmetry=crystal_symmetry,
         dummy_refinement=dummy_refinement,
         refinement_params=params.input_files.refinement_params,
         correct_special_position_tolerance=
           params.non_user_params.correct_special_position_tolerance,
         ncs=ncs,
         skip_clash_guard=params.non_user_params.skip_clash_guard ,
         remove_clashing_residues=\
           params.refine_top_models.remove_clashing_residues,
         params=params,Facts={},
         clash_cutoff=params.refine_top_models.clash_cutoff) # 2013-09-12
        self.wait_for_file_to_appear(refined_mr_model,
         max_ticks=params.control.max_wait_time,allow_failure=False)

        print("Copying data from %s to %s and setting space_group= %s" %(
         params.input_files.data,new_data,
         self.get_space_group_symbol(crystal_symmetry=crystal_symmetry)), file=out)
        self.copy_to_refine_data(params.input_files.data,new_data,
         crystal_symmetry=crystal_symmetry,out=out)

      if params.place_model.denmod_after_refine and \
           params.place_model.refine_after_mr:
        print("\nGetting density modified map coefficients...", file=out)
        map_coeffs=self.build_with_autobuild(params,solution=None,
          out=out,data=new_data,maps_only=True,
          start_map_coeffs=map_coeffs,
          maps_only_model=refined_mr_model,
          ps_in_rebuild=params.place_model.ps_in_rebuild)
        labin_map_coeffs="FP=FWT PHIB=PHWT"

      refined_mr_model_list.append(refined_mr_model)
      refined_mr_model_mtz_list.append(map_coeffs)
      data_list.append(new_data)
    return refined_mr_model_list,refined_mr_model_mtz_list,\
         data_list,crystal_symmetry_list,labin_map_coeffs

  def copy_to_refine_data(self,data,output_data,crystal_symmetry=None,
     out=sys.stdout):
    from phenix.autosol.merge_mtz import merge_mtz
    arrays_list=[["ALL"]] # take everything and set crystal symmetry
    file_list=[data]
    merge_mtz=merge_mtz(file_list=[data],arrays_list=[["ALL"]],
       output_file=output_data,crystal_symmetry=crystal_symmetry,out=out,
       skip_phases=True) # don't copy PHIM because it is not mapped with new sym
    merge_mtz.select_arrays()
    merge_mtz.write_arrays()

  def run_comparison(self,params,solution=None,out=sys.stdout):
    model=solution.model
    map_coeffs=solution.map_coeffs
    labin_map_coeffs=solution.labin_map_coeffs
    print("\n",80*"-", file=out)
    print("\nComparison of model/map with comparison file: ", file=out)
    if model is not None:
      print("\nModel:%s " %(model), file=out)
    if model is not None:
      map_cc,map_cc_local=self.get_model_comparison_cc(
           params,model=model,out=out)
      if map_cc is not None:
        solution.model_cc_to_comparison=map_cc
        solution.model_cc_local_to_comparison=map_cc_local
    if map_coeffs is not None:
      print("\nMap:%s  labin:%s" %(map_coeffs,labin_map_coeffs), file=out)
      cc=self.get_map_comparison_cc(
          params,map_coeffs=map_coeffs,
          labin_map_coeffs=labin_map_coeffs,out=out)
      if cc is not None:
        solution.cc_to_comparison=cc
    print("\n",80*"-", file=out)

  def get_model_comparison_cc(self,params,model=None,out=sys.stdout):
    mtz=params.non_user_params.comparison_mtz
    labin_mtz=params.non_user_params.labin_comparison_mtz
    from phenix.command_line.get_cc_mtz_pdb import get_cc_mtz_pdb
    quiet=not params.control.verbose
    from six.moves import cStringIO as StringIO
    f=StringIO()
    args=["mtz_in="+str(mtz),
        "pdb_in="+str(model),
        "quick=True",
        "temp_dir="+str(params.directories.temp_dir),
        "output_dir="+str(params.directories.temp_dir)]
    if labin_mtz: args.append("labin="+labin_mtz)
    g=get_cc_mtz_pdb(args,quiet=quiet,out=f)
    if g.found_overall is None or g.found_region is None:
      print("No model correlation obtained", file=out)
      return None,None
    else:
      print("Model correlation is %6.3f (overall) " %(g.found_overall) +\
       "and %6.3f (in region of model)" %(g.found_region), file=out)
      return g.found_overall,g.found_region

  def get_map_comparison_cc(self,params,map_coeffs=None,labin_map_coeffs=None,
       out=sys.stdout):
    mtz=params.non_user_params.comparison_mtz
    labin_mtz=params.non_user_params.labin_comparison_mtz
    cc=self.get_map_cc(params,map1=mtz,labin_1=labin_mtz,
      map2=map_coeffs,labin_2=labin_map_coeffs,out=sys.stdout)
    if cc is None:
      print("No map correlation obtained ", file=out)
      return None
    else:
      print("Map correlation is %6.3f" %(cc), file=out)
      return cc

  def remove_free(self,mtz=None,free=None,new_mtz=None,out=sys.stdout,
     labin="FP=FWT PHIB=PHWT",labin_free="FreeR_flag=FreeR_flag",
     temp_dir="",test_flag_value=None):
    print("\nRemoving free reflection data from %s to yield %s " %(
       mtz,new_mtz), file=out)
    args=["map_coeffs="+str(mtz),
          "free_in="+free,
          "mtz_out="+new_mtz,
          "temp_dir="+temp_dir,
         ]
    if test_flag_value:
       args.append("test_flag_value="+test_flag_value)
    if labin_free:
       args.append("labin_free="+labin_free)
    if labin:
       args.append("labin_map_coeffs="+labin)

    log=os.path.join(new_mtz[:-4]+".log")
    f=open(log,'w')
    from phenix.command_line.remove_free_from_map import remove_free_from_map
    remove_free_from_map(args,out=f)
    f.close()
    new_labin="FP=FWT PHIB=PHWT"
    print("Log file is %s" %(log), file=out)
    return new_labin

  def convert_to_map(self,params,mr_model_mtz_list,labin="FP=FWT PHIB=PHWT",
       out=sys.stdout,remove_free=True):
    #NOTE: writes to same directory as mtz files are in
    mr_model_map_list=[]
    for mtz in mr_model_mtz_list:
      if mtz is None:
        mr_model_map_list.append(None)
      else:
        from phenix.autosol.trim_file_name import trim_file_name
        if (not params.control.real_space_optimize) and remove_free and \
           (not hasattr(params.input_files,'remove_free') or \
                    params.input_files.remove_free):
          new_mtz=mtz[:-4]+"_nf.mtz"
          new_labin=self.remove_free(mtz=mtz,free=params.input_files.data,
            labin_free=params.input_files.labin,
            new_mtz=new_mtz,labin=labin,
            test_flag_value=params.control.test_flag_value,
            out=out,temp_dir=params.directories.temp_dir)
        else:
          new_mtz=mtz
          new_labin=labin
          if new_labin is None:
            raise Sorry("Need labin for map coeffs defined if map_coeffs are supplied and free reflections are not removed. Use labin_map_coeffs='FP=FWT PHIB=PHWT' or equivalent")
        self.wait_for_file_to_appear(new_mtz,
           max_ticks=params.control.max_wait_time,allow_failure=False)
        new_mtz_trim=trim_file_name(new_mtz).trimmed_file
        full_new_mtz_trim=os.path.join(params.directories.temp_dir,new_mtz_trim)
        self.copyfile(new_mtz,full_new_mtz_trim)
        self.wait_for_file_to_appear(full_new_mtz_trim,
           max_ticks=params.control.max_wait_time,allow_failure=False)
        map=new_mtz[:-4]+".map"
        map_trim=new_mtz_trim[:-4]+".map"

        text="ccp4_map_file "+str(map_trim)+"\n"
        text+="mask_cycles 1 \n minor_cycles 0\nsolvent_content 0.5\n no_build\n"
        from phenix.autosol.run_resolve import run_resolve
        resolve=run_resolve(
          temp_dir=params.directories.temp_dir,
          hklin=new_mtz_trim,
          labin=new_labin,
          build="no_build",
          command_1=text)
        if not os.path.isfile(os.path.join(params.directories.temp_dir,map_trim)):
          raise Sorry("No map file produced for %s" %(map))
        print("Converted %s to map using labels %s " %(
             mtz,labin), file=out)
        full_map_trim=os.path.join(params.directories.temp_dir,map_trim)
        self.wait_for_file_to_appear(full_map_trim,
           max_ticks=params.control.max_wait_time,allow_failure=False)
        self.copyfile(full_map_trim,map)
        self.wait_for_file_to_appear(map,
           max_ticks=params.control.max_wait_time,allow_failure=False)
        mr_model_map_list.append(map)
    return mr_model_map_list

  def score_models(self,params=None,models=[],out=sys.stdout,
     mtz_in=None,labin=None,temp_dir="",write_map_coeffs=False,
      apply_ncs=True,ncs_object=None,crystal_symmetry=None):
    # allow crystal_symmetry explicitly  if always same
    # if write_map_coeffs=true and only one model, write and return map coeffs file
    print("\nScoring %d models" %(len(models)), file=out)
    if write_map_coeffs:
       print("Writing map coeffs", file=out)
    if temp_dir:
      print("Working directory for scoring: %s \n" %(temp_dir), file=out)

    self.results=[]
    best_llg_score=None
    id=0
    llg_score_list=[]
    for model in models:
      print("Scoring %s" %(model), file=out)
      id+=1
      if mtz_in is None: mtz_in=params.input_files.mtz_in
      if labin is None: labin=params.input_files.labin
      self.wait_for_file_to_appear(model,
       max_ticks=params.control.max_wait_time,allow_failure=False)
      self.wait_for_file_to_appear(mtz_in,
       max_ticks=params.control.max_wait_time,allow_failure=False)

      labin_fp_sigfp=self.get_labels_from_labin_string(
         labin,selection=['FP','SIGFP']).split()
      if len(labin_fp_sigfp)<2:
        raise Sorry("Labin line needs FP and SIGFP for mr_rescoring "+str(labin))
      identity=params.place_model.identity_for_scoring_only
      # standard value not actual

      for f,name in zip(
        [model,mtz_in],
        ['model','mtz_in']):
        if not os.path.exists(f) :
          raise Sorry("\n%s (%s) seems to be missing?" %(f,name))

      # explicitly specify crystal symmetry from model
      pdb_input = get_pdb_inp(file_name=model)
      if crystal_symmetry is None:
        if hasattr(self,'crystal_symmetry') and \
             self.crystal_symmetry is not None:
          crystal_symmetry=self.crystal_symmetry
        else:
          crystal_symmetry=self.none_if_no_space_group_info(
            pdb_input.crystal_symmetry_from_cryst1())
      if crystal_symmetry is None:
        raise Sorry("No crystal symmetry found in %s " %(model))

      if params.rescore_mr.edit_model:
        # here if we come in with a rosetta model that needs generation of entire pdb

        model=self.apply_model_info(params,model=model,
          apply_ncs=apply_ncs,ncs_object=ncs_object,out=out,
          crystal_symmetry=crystal_symmetry)

      if params.non_user_params.file_base:
        root=os.path.join(temp_dir,params.non_user_params.file_base)
      else:
        root=self.get_file_base(model)
      if write_map_coeffs:
        map_coeffs=root+".1.mtz"
      else:
        map_coeffs=None
      score_log=os.path.join(params.directories.temp_dir,
          "%s_score_%d.log" % (params.non_user_params.file_base,id))
      from phaser import InputMR_DAT,runMR_DAT,InputMR_RNP,runMR_RNP
      if 1:
        i = InputMR_DAT()
        if params.control.verbose:
          i.setMUTE(False)
        else:
          i.setMUTE(True)
        i.setHKLI(mtz_in)
        i.setLABI_F_SIGF(labin_fp_sigfp[0],labin_fp_sigfp[1])
        r = runMR_DAT(i)
        if not r.Success():
           if not params.control.verbose:
             print("Phaser rescoring failed... rerunning to show log", file=out)
             i.setMUTE(False)
             r = runMR_DAT(i)
             print("End of failed phaser rescoring", file=out)
           raise Sorry("\nSorry phaser rescoring failed with the model %s " %(
             model) )
        i = InputMR_RNP()
        i.setJOBS(1) # if compiled with OpenMP assume load balancing is done elsewhere
        if params.control.verbose:
          i.setMUTE(False)
        else:
          i.setMUTE(True)
        if hasattr(r,'getTopSet'):
          i.setSPAC_HALL(r.getTopSet().getSpaceGroupHall)
        else:
          i.setSPAC_HALL(r.getSpaceGroupHall())
        i.setCELL6(r.getUnitCell())
        i.setREFL_F_SIGF(r.getMiller(),r.getFobs(),r.getSigFobs())
        i.setROOT(root)
        i.addENSE_PDB_ID("model",model,float(identity)) # XXX can we add model as text?
        i.addSOLU_ORIG_ENSE("model")
        if write_map_coeffs:
          i.setHKLO(True)
        else:
          i.setHKLO(False)
        r = runMR_RNP(i)
        if not r.Success() or not r.foundSolutions():
           if not params.control.verbose:
             print("Phaser rescoring failed...rerunning to show log", file=out)
             i.setMUTE(False)
             r = runMR_RNP(i)
             print("End of failed phaser rescoring", file=out)
           raise Sorry("\nSorry phaser rescoring failed with the model %s " %(
             model) )
        llg=r.getTopLLG()
      else:
        llg=10.
        print("SKIPPING MR SCORING")
      llg_score_list.append(llg)
      if best_llg_score is None or llg > best_llg_score: best_llg_score=llg
      print("Score for %s is %6.2f \n"% (model,llg), file=out)
      if map_coeffs:
         print("Map coefficients file: %s " % (map_coeffs), file=out)
         labin_map_coeffs="FP=FWT PHIB=PHWT"
      else:
         labin_map_coeffs="FP=FWT PHIB=PHWT"
      from phenix.rosetta.mr_rosetta import mr_rosetta_solution
      if hasattr(self,'highest_id'):
        self.highest_id+=1
        id=self.highest_id
      else:
        id=1
        self.highest_id=id
      solution=mr_rosetta_solution(name=model,model=model,id=id,
          stage="scored_solution",mr_llg=llg,map_coeffs=map_coeffs,
          labin_map_coeffs=labin_map_coeffs,crystal_symmetry=crystal_symmetry)
      self.results.append(solution)
    print("LLG scores: ", end=' ', file=out)
    for llg in llg_score_list:
      print(" %8.2f " %(llg), end=' ', file=out)
    if best_llg_score is not None:
      print("\nDone with scoring this model . Best score =%6.2f" %(
        best_llg_score), file=out)
    else:
      print("\nDone with scoring this model. No LLG score obtained", file=out)
    return best_llg_score

  def setpaths(self,params,out=sys.stdout,check_paths=False):
    if params.control.verbose and hasattr(self,'rosetta_script_dir'):
      print("Setting paths", file=out)
      print("Rosetta scripts are in: ",self.rosetta_script_dir, file=out)
      print("Rosetta is in: ",self.rosetta_path, file=out)
    if check_paths:
      if not os.path.isdir(self.rosetta_script_dir):
        raise Sorry("Rosetta script dir is missing %s" \
            %(str(self.rosetta_script_dir)))
      if not os.path.isfile(self.rosetta_path):
        raise Sorry("Rosetta is missing %s" %(str(self.rosetta_path)))

  def get_solution_key(self,input_pdb,allow_fewer=False):
    # take 1st 5 chars of file name (2010-10-27 was 4, now 5 to match rosetta)
    from phenix.autosol.trim_file_name import trim_file_name
    text=trim_file_name(input_pdb).trimmed_file
    if not allow_fewer and len(text[:-4])<4:
      raise Sorry("Sorry, input PDB file names for mr_rosetta need to be 4 chars plus '.pdb'")
    return text[:5]

  def divisors(self,n): # return divisors of n
    dd=[]
    for i in range(1,n+1):
      j=n//i
      if j*i == n:
        dd.append(i)
    return dd

  def read_hhr_files(self,params,out=sys.stdout): # skip if alignment files
    hhr_file_list=params.input_files.hhr_files
    if params.input_files.alignment_files or not hhr_file_list or \
       params.place_model.model_already_aligned: return

    if params.place_model.model_already_placed:  # not allowed
      raise Sorry(
       "Sorry, you cannot use model_already_placed=True with hhr_files")

    # figure out how many copies of model we want to try to extract
    if not params.read_hhpred.copies_to_extract:
      print("\nGuessing range of copies_to_extract using ncs_copies=%s" %(
       str(params.crystal_info.ncs_copies)), file=out)
      copies_to_extract=[1]
      for ncs_copies in params.crystal_info.ncs_copies:
        if ncs_copies is None: continue
        for value in self.divisors(ncs_copies):
          if not value in copies_to_extract:
            copies_to_extract.append(value)
      copies_to_extract.sort()
      params.read_hhpred.copies_to_extract=copies_to_extract
    print("Copies to extract: %s" %(
      str(params.read_hhpred.copies_to_extract)), file=out)

    params.input_files.search_models=[]
    params.input_files.copies_in_search_models=[]
    params.input_files.alignment_files=[]

    for hhr_file in hhr_file_list:
      from phenix.autosol.trim_file_name import trim_file_name
      text=trim_file_name(hhr_file).trimmed_file[:-4]
      workdir=self.create_temp_directory(n_dir_max=params.control.n_dir_max,
        path=params.directories.temp_dir,prefix=text)
      print("Working directory for set of alignment files will be: %s " %(workdir), file=out)
      edited_hhr_file=os.path.join(workdir,text+"_ed.hhr")
      entry_name_list=self.edit_hhr_file(params,hhr_file=hhr_file,
         edited_hhr_file=edited_hhr_file,out=out)
      edited_pdb_file_list,alignment_file_list,copies_list=self.read_hhr_file(
        params,hhr_file=edited_hhr_file,workdir=workdir,
        entry_name_list=entry_name_list,
        copies_to_extract_list=params.read_hhpred.copies_to_extract,out=out)

      # now append these edited PDB and alignment files...
      params.input_files.search_models+=edited_pdb_file_list
      params.input_files.copies_in_search_models+=copies_list
      params.input_files.alignment_files+=alignment_file_list
    params.place_model.model_already_aligned=True

  def parse_hhr_file(self,hhr_file=None,out=sys.stdout):
    # split hhr file into pieces, one for each pdb
    header=""
    data_text=""
    entry_info_list=[]
    entry_name_list=[]
    entry_data_list=[]
    finished_header=False
    finished_info_list=False

    for line in open(hhr_file).readlines():
      line=line.rstrip()
      if line.find("Done")==0:
        break # done
      if line.find(" No Hit")==0: # start of info lines
        header+=line+"\n"
        if not finished_header:
          finished_header=True

      elif finished_header and line.find("No ")==0: # start or continue data lines
        if not finished_info_list:
          finished_info_list=True
        if data_text:
          entry_data_list.append(data_text.rstrip())
          data_text=""
        data_text+=line+"\n"

      elif not finished_header:
        header+=line+"\n"
      elif not finished_info_list:
        if line:
          entry_info_list.append(line.rstrip())
        else: pass # skip blank lines in info list
      else: #continue on current data list
        data_text+=line+"\n"
    if data_text:
      entry_data_list.append(data_text.rstrip())

    entry_name_list=[]
    for line in entry_info_list:
      try:
        pdb_name=line.split()[1][:4] # looks like "32 1g7n_A Adipocyte lipid_binding..."
      except Exception:
        raise Sorry("Sorry cannot interpret the line \n%s\n from %s" %(line,hhr_file))
      entry_name_list.append(pdb_name)

    if len(entry_data_list) < len(entry_info_list):
      print("NOTE: Only using the first %d entries in hhr file %s" %(
      len(entry_data_list),hhr_file), file=out)
      entry_info_list=entry_info_list[:len(entry_data_list)]
    elif len(entry_data_list) > len(entry_info_list):
      raise Sorry("Sorry, number of pdb entries in %s (%d) does not match " %(
        hhr_file,len(entry_data_list)) +\
         "\nnumber of alignnment entries (%d)" %(len(entry_info_list)))
    return header,entry_info_list,entry_data_list,entry_name_list

  def edit_hhr_file(self,params,hhr_file=None,edited_hhr_file=None,
       workdir=None,out=sys.stdout):
    number_of_models_to_skip=params.read_hhpred.number_of_models_to_skip
    number_of_models=params.read_hhpred.number_of_models
    print("\nEditing the hhr file %s -> %s " %(hhr_file,edited_hhr_file), file=out)
    if number_of_models:
      print("Taking up to the first %d models from hhr file " %(
        number_of_models), file=out)
    if number_of_models_to_skip:
      print("(After first skipping first %d models)" %(
         number_of_models_to_skip), file=out)
    # edit the file to take out entries we don't want
    header,entry_info_list,entry_data_list,entry_name_list=self.parse_hhr_file(
        hhr_file=hhr_file,out=out)
    edited_pdb_list=[]
    pdb_name_list=[]

    if number_of_models_to_skip:
      entry_info_list=entry_info_list[number_of_models_to_skip:]
      entry_data_list=entry_data_list[number_of_models_to_skip:]
      entry_name_list=entry_name_list[number_of_models_to_skip:]
    if number_of_models:
      entry_info_list=entry_info_list[:number_of_models]
      entry_data_list=entry_data_list[:number_of_models]
      entry_name_list=entry_name_list[:number_of_models]

    f=open(edited_hhr_file,'w')
    print(header.rstrip(), file=f)
    for entry in entry_info_list:
      print(entry, file=f)
    print(file=f) # blank line
    for entry in entry_data_list:
      print(entry+"\n\n", file=f)
    print("Done!\n", file=f)
    f.close()
    self.wait_for_file_to_appear(edited_hhr_file,
       max_ticks=params.control.max_wait_time,allow_failure=False)
    return entry_name_list

  def read_hhr_file(self,params,hhr_file=None,workdir=None,entry_name_list=None,
        copies_to_extract_list=[1], out=sys.stdout):
    print("\nLoading PDB and alignment files from %s" %(hhr_file), file=out)
    if params.place_model.align_with_sculptor:
       edited_pdb_file_list,alignment_file_list,copies_list=\
          self.read_hhr_file_with_phenix(params,hhr_file=hhr_file,
        workdir=workdir,entry_name_list=entry_name_list,
        copies_to_extract_list=copies_to_extract_list, out=out)
    else:
       print("\nUsing rosetta scripts to read hhr files "+\
         "(mr_model_preparation disabled)\n", file=out)
       edited_pdb_file_list,alignment_file_list,copies_list=\
          self.read_hhr_file_with_rosetta(params,hhr_file=hhr_file,
        workdir=workdir,entry_name_list=entry_name_list,
        copies_to_extract_list=copies_to_extract_list, out=out)
    return edited_pdb_file_list,alignment_file_list,copies_list

  def read_hhr_file_with_phenix(self,params,hhr_file=None,workdir=None,
     entry_name_list=None,copies_to_extract_list=[1], out=sys.stdout):
    # NOTE: this aligns the file, so no alignment file is necessary
    log_file=os.path.join(workdir,'mr_model_preparation.log')
    args=["file_name="+str(hhr_file),
          "max_hits="+str(params.read_hhpred.number_of_models),
          "remove_alternate_conformations=True",
          "rename=True",
          "output.folder="+str(workdir)]
    f=open(log_file,'w')
    print("Log file for mr_model_preparation is: %s "%(log_file), file=out)
    from phaser import mr_model_preparation
    file_list=mr_model_preparation.run(args=args,out=f)
    f.close()
    edited_pdb_file_list=[]
    alignment_file_list=[]
    copies_list=[]
    for entry,file_names in zip(entry_name_list,file_list):
      for copies_to_extract in copies_to_extract_list:
        if len(file_names)<copies_to_extract:
          continue
        print("Extracting %d copies..."%(copies_to_extract), file=out)
        text=""
        for file_name in file_names[:copies_to_extract]:
          text+=open(file_name).read()
        edited_pdb_file=os.path.join(workdir,entry+"_"+str(copies_to_extract)+
           "_chain.pdb")
        f=open(edited_pdb_file,'w')
        print(text, file=f)
        f.close()
        alignment_file='None'
        self.wait_for_file_to_appear(edited_pdb_file)

        edited_pdb_file_list.append(edited_pdb_file)
        copies_list.append(copies_to_extract)
        alignment_file_list.append(None)
    print("Search model list: %s\n" %(edited_pdb_file_list), file=out)
    return edited_pdb_file_list,alignment_file_list,copies_list

  def read_hhr_file_with_rosetta(self,params,hhr_file=None,workdir=None,
     entry_name_list=None,copies_to_extract_list=[1], out=sys.stdout):
    rosetta_script=os.path.join(
       self.rosetta_script_dir,"prepare_template_for_MR.pl")
    if not os.path.isfile(rosetta_script):
      raise Sorry("\nThe script file %s seems to be missing?" %(rosetta_script))

    run_file=os.path.join(workdir,'read_hhr_file.sh')
    log_file=os.path.join(workdir,'read_hhr_file.log')

    commands="#!/bin/sh  \n"
    commands+="cd %s\n" %(self.add_double_quote(workdir,escape_space=False))
    commands+="%s %s > %s \n" % (
      self.add_double_quote(rosetta_script,escape_space=False),
      self.add_double_quote(hhr_file),self.add_double_quote(log_file,escape_space=False))
    if params.control.verbose:
      print("Commands for alignment: \n %s" %(commands), file=out)
    f=open(run_file,'w')
    print(commands, file=f)
    f.close()
    os.chmod(run_file,0o755)
    try:
      self.start_run(run_file)
    except Exception as e:
      if str(e).find("ftp.wwpdb.org"):  # assume we could not connect
        line="\n"+str(e)+\
         "\nSorry, could not connect to ftp.wwpdb.org to download PDB files..."+\
         "\nIf you already have your edited pdb files (e.g., 1qor_mr.pdb) and "+\
         "alignment files\n(e.g., 1qor.ali) you can just enter them with "+\
         "alignment_files=1qor.ali  \nsearch_models=1qor_mr.pdb \n"
        raise Sorry(line)

    edited_pdb_file_list=[]
    alignment_file_list=[]
    copies_list=[]
    copies=1 # always
    for entry in entry_name_list:
      edited_pdb_file=os.path.join(workdir,entry+"_mr.pdb")
      alignment_file=os.path.join(workdir,entry+".ali")
      self.wait_for_file_to_appear(edited_pdb_file)
      self.wait_for_file_to_appear(alignment_file)
      further_edited_pdb_file=os.path.join(workdir,entry+"_mr_ed.pdb")
      final_file=self.remove_het(edited_pdb_file,
           further_edited_pdb_file,out=out)
      edited_pdb_file_list.append(final_file)
      alignment_file_list.append(alignment_file)
    print("EDITED PDB FILE LIST: %s" %(edited_pdb_file_list), file=out)
    print("ALIGNMENT FILE LIST: %s" %(alignment_file_list), file=out)
    copies_list.append(copies)
    return edited_pdb_file_list,alignment_file_list,copies_list

  def remove_het(self,input_file,output_file,out=sys.stdout):
    # just remove HETATM
    f=open(output_file,'w')
    pdb_input = get_pdb_inp(file_name=input_file)
    crystal_symmetry=self.none_if_no_space_group_info(
       pdb_inp.crystal_symmetry_from_cryst1())
    if crystal_symmetry is None and hasattr(self,'crystal_symmetry'):
      crystal_symmetry=self.crystal_symmetry
    if crystal_symmetry is not None:
      print(iotbx.pdb.format_cryst1_record(
           crystal_symmetry=crystal_symmetry), file=f)

    pdb_hierarchy=pdb_inp.construct_hierarchy()
    has_het=False
    for model in pdb_hierarchy.models()[:1]:
      for chain in model.chains():
        for conformer in chain.conformers()[:1]:
          for residue in conformer.residues():
            atoms=residue.atoms()
            for atom in atoms:
              if atom.hetero:
                has_het=True
              else:
                 print(atom.format_atom_record(), file=f)
    f.close()
    if has_het:
      print("\nCopying %s to %s, removing HETATM records\n" %(
       input_file,output_file), file=out)
      return output_file
    else:
      return input_file

  def simple_write_one_alignment_file(self,chain_id=None,
        out=sys.stdout,base_name="coords"):
      # return text for one of the component alignment files
      group=self.group_by_chain_id[chain_id]
      seq_file_start=self.seq_file_start_by_group[group]
      pdb_file_start=self.pdb_file_start_by_group[group]
      pdb_file_alignment=self.pdb_file_alignment_by_group[group]
      seq_file_alignment=self.seq_file_alignment_by_group[group]
      from six.moves import cStringIO as StringIO
      f=StringIO()
      print("## TARGET %s " %(base_name), file=f)
      print("# hhsearch", file=f)
      print("scores_from_program: 0 1.00", file=f)
      print("%d %s" %(seq_file_start,seq_file_alignment), file=f)
      print("%d %s" %(pdb_file_start,pdb_file_alignment), file=f)
      print("--", file=f)
      return f.getvalue()

  def has_non_het(self,chain):
        has_non_het=False
        for conformer in chain.conformers()[:1]:
          for residue in conformer.residues():
            atoms=residue.atoms()
            for atom in atoms:
              if not atom.hetero:
                has_non_het=True
        return has_non_het

  def split_input_pdb_by_group(self,params,input_pdb=None,out=sys.stdout):
    # writes model into separate PDB files for each chain type (group) and
    # sets ncs_object_by_group,model_by_group

    # Use to split up model for running rosetta on separate chains.
    # Follow with rosetta, then apply_model_info using
    # ncs_object_by_group to apply ncs for this set, then combine resulting files
    # to yield all chains and all ncs copies 2010-12-11

    self.read_model(input_pdb,params=params,out=out,quiet=True)
    #  so we can have group_from_chain_id

    pdb_input = get_pdb_inp(file_name=input_pdb)
    hierarchy = pdb_input.construct_hierarchy()

    # requires self.group_by_chain_id self.seq_file_by_group self.align_file_by_group,
    ncs_object_by_group=[]
    model_by_group=[]
    group_by_chain_id={}
    file_list=[]
    have_something={}

    base_file=os.path.join(params.directories.temp_dir,
      os.path.split(input_pdb)[-1][:-4])
    for group in range(self.number_of_groups):
      model=base_file+"_"+str(group)+".pdb"
      model_by_group.append(model)
      file=open(model,'w')
      file_list.append(file)
      print(self.cryst1_record_from_anywhere(pdb_input=pdb_input), file=file)
      print("Writing chains for group %s to %s" %(str(group),file.name), file=out)

    for model in hierarchy.models()[:1]:
      for chain in model.chains():
        if not self.has_non_het(chain): continue # skip this chain
        chain_id=chain.id.replace(" ","")
        if not chain_id: chain_id="A"
        group=self.get_group_from_chain_id(chain_id,out=out)
        group_by_chain_id[chain_id]=group
        for conformer in chain.conformers()[:1]:
          for residue in conformer.residues():
            atoms=residue.atoms()
            for atom in atoms:
              if atom.hetero:
                has_het=True
              else:
                print(atom.format_atom_record(), file=file_list[group])
                have_something[group]=True

    for group in range(self.number_of_groups):
      file_list[group].close()
      if not group in have_something.keys() or not have_something[group]:
        raise Sorry("\nSorry, the file %s representing group %d has no atoms? " %(
            file_list[group].name,group)+
            "\nThis could indicate a problem with sequence alignment in rebuilding.\n"+
            "It could be that there are near-duplicate sequences in your sequence file\n"+
            "...if this is the case then you may need to artificially rename your duplicated\n"+
            "chains and renumber residues to make a single chain containing a duplication\n" )
      self.wait_for_file_to_appear(model_by_group[group],
          max_ticks=params.control.max_wait_time,allow_failure=False)
      # get an ncs object for each group now
      ncs_object=self.get_ncs_from_mr_models(params,
        [model_by_group[group]],out=out)[0] # goes in and comes back as a list
      ncs_object_by_group.append(ncs_object)
    return ncs_object_by_group,model_by_group,group_by_chain_id

  def simple_split_input_pdb(self,params,out=sys.stdout):
    hierarchy = self.pdb_input.construct_hierarchy()
    print("\nSplitting model into files for each chain\n", file=out)
    model_text_list=[]
    chain_id_list=[]
    for model in hierarchy.models()[:1]:
      for chain in model.chains():
        # determine if chain has any non-het atoms
        if not self.has_non_het(chain): continue # skip this chain
        chain_id=chain.id.replace(" ","")
        if not chain_id: chain_id="A"
        from six.moves import cStringIO as StringIO
        f=StringIO()
        print(self.cryst1_record_from_anywhere(pdb_input=self.pdb_input), file=f)
        for residue_group in chain.residue_groups():
          for atom_group in residue_group.atom_groups():
            atom_group=self.replace_mse_in_atom_group(atom_group)
            for atom in atom_group.atoms():
              print(atom.format_atom_record(), file=f)

        print("Converted chain %s " %(chain_id), file=out)
        text=f.getvalue()
        chain_id_list.append(chain_id)
        model_text_list.append(text)
    return model_text_list,chain_id_list

  def replace_mse_in_atom_group(self,atom_group):
    if atom_group.resname.lower().replace(" ","")=='mse':
      atom_group.resname='MET'
      for atom in atom_group.atoms():
        atom.hetero=False
        if atom.name.lower().replace(" ","")=='se':
          atom.name=" SD"
          atom.element="S"
    return atom_group

  def cryst1_record_from_anywhere(self,pdb_input=None):
        crystal_symmetry=None
        if pdb_input:
          crystal_symmetry=self.none_if_no_space_group_info(
            pdb_input.crystal_symmetry_from_cryst1())
        if not crystal_symmetry and \
            hasattr(self,'crystal_symmetry') and self.crystal_symmetry:
          crystal_symmetry=self.crystal_symmetry
        if not crystal_symmetry:
          crystal_symmetry=self.dummy_crystal_symmetry()
        if not hasattr(self,'crystal_symmetry') or not self.crystal_symmetry:
          self.crystal_symmetry=crystal_symmetry # save it
        return iotbx.pdb.format_cryst1_record(
          crystal_symmetry=crystal_symmetry)+"\n"

  def simple_create_aligned_model(self,params,out=sys.stdout):
    # apply alignment to a single PDB file.
    # also create new simple alignment file and dict of B values
    if hasattr(self,'model_info'): return

    input_pdb=params.input_files.search_models[0]

    if not hasattr(params,'place_model'):
      already_aligned=True
      alignment_file=None
    elif params.place_model.model_already_aligned:
      already_aligned=True
      alignment_file=None
      print("\nMODEL IS ALREADY ALIGNED TO SEQUENCE FILE\n", file=out)
    elif params.input_files.alignment_files:
       alignment_file=params.input_files.alignment_files[0]
       already_aligned=False
    else:
       print("No alignment file supplied...using search model as is", file=out)
       already_aligned=True
       alignment_file=None
       print("\nMODEL IS NOT ALIGNED TO SEQUENCE FILE\n", file=out)

    # make sure we have 5-char name:
    input_pdb=self.check_5_chars_pdb(params,input_pdb,alignment_file)

    try_to_get_symmetry=self.decide_to_get_symmetry(params) # True if placed

    self.read_seq_file(params,out=out) # now self.sequence=[seq1,seq2...]
    if params.control.verbose:
       print("\nSequence(s) from sequence file: ",self.sequence, file=out)
    self.simple_read_alignment_file_and_model(params,
      alignment_file=alignment_file,model=input_pdb,out=out,
       try_to_get_symmetry=try_to_get_symmetry)
    self.solution_key=self.get_solution_key(input_pdb)
    params.non_user_params.file_base=self.solution_key
    if already_aligned or (
        params.place_model.force_alignment and params.input_files.alignment_files):
      params.place_model.model_already_aligned=True
      return input_pdb

    # Apply alignment to input_pdb file

    # 2010-12-09 we want to go through each chain in input_pdb, decide which
    # group it is from its chain_id, apply the appropriate alignment file
    # then assemble them all at the end

    aligned_search_model=os.path.join(params.directories.temp_dir,
       os.path.split(input_pdb)[-1]+"_al.pdb")
    print("\nApplying alignment in %s to input search model %s to \nyield %s\n" \
      %(str(alignment_file),input_pdb,aligned_search_model), file=out)

    ncs_copies=params.crystal_info.ncs_copies[0]

    # split up model, apply alignments, recombine to make new aligned model

    model_text_list,chain_id_list=self.simple_split_input_pdb(params,
        out=out)
    text=self.cryst1_record_from_anywhere(pdb_input=self.pdb_input)

    for model_text,chain_id in zip(model_text_list,chain_id_list):
      alignment_text=self.simple_write_one_alignment_file(
        chain_id=chain_id,out=out)
      print("ALIGNMENT TEXT: \n",alignment_text)
      text+=self.simple_align_with_sculptor(params,chain_id=chain_id,
          alignment_text=alignment_text,model_text=model_text,out=out)
    alignment_file=os.path.join(params.directories.temp_dir,
        'alignment.ali')
    f=open(alignment_file,'w')
    print(text, file=f)
    f.close()
    self.wait_for_file_to_appear(alignment_file,
        max_ticks=params.control.max_wait_time,allow_failure=False)
    print("\nAlignment written to %s " %(alignment_file), file=out)
    params.input_files.alignment_files=[alignment_file]
    params.place_model.model_already_aligned=True

    aligned_search_model=os.path.join(params.directories.temp_dir,
      os.path.split(input_pdb)[-1][:-4]+"_al.pdb")
    f=open(aligned_search_model,'w')
    print(text, file=f)
    f.close()
    self.wait_for_file_to_appear(aligned_search_model,
       max_ticks=params.control.max_wait_time,allow_failure=False)
    print("Wrote aligned search model to %s" %(
      aligned_search_model), file=out)


    # and now finally redo everything with this aligned model
    print("\nRerunning PDB read with aligned model\n", file=out)
    params.input_files.search_models=[aligned_search_model]
    delattr(self,'model_info')
    aligned_search_model=\
     self.simple_create_aligned_model(params,out=sys.stdout) # to get model info
    return aligned_search_model

  def check_5_chars_pdb(self,params,search_model,alignment_file):
    from phenix.autosol.trim_file_name import trim_file_name
    text=trim_file_name(search_model).trimmed_file[:-4]
    if len(text)<5:
      if alignment_file: # can't do it
        raise Sorry("Sorry, if an alignment file is provided then all "+
         "search models \nmust have at least 5 characters before the '.pdb'" )
      else:  # copy the file to one with more characters
        new_file=os.path.join(self.params.directories.workdir,text+"_copy.pdb")
        new_file=self.make_full_path(new_file)
        self.copyfile(search_model,new_file)
        new_file=self.make_full_path(new_file)
      return new_file
    else: # it's fine already
      return search_model

  def simple_align_with_sculptor(self,params,alignment_text="",
         out=sys.stdout,chain_id=None,model_text="",
    pdb_file_alignment=None,
    seq_file_alignment=None,):
    print("Aligning with sculptor", file=out)

    # edit alignment file to be in format for sculptor: 2010-11-30
    if pdb_file_alignment is None or seq_file_alignment is None:
      group=self.group_by_chain_id[chain_id]
      if pdb_file_alignment is None:
        pdb_file_alignment=self.pdb_file_alignment_by_group[group]
      if seq_file_alignment is None:
        seq_file_alignment=self.seq_file_alignment_by_group[group]

    align_text=">sequence\n%s\n>pdb\n%s" %(
        seq_file_alignment,pdb_file_alignment)

    #JUST XXX pass as strings if possible
    pdb_file=os.path.join(params.directories.temp_dir,"sculpt.pdb")
    f=open(pdb_file,'w')
    print(model_text, file=f)
    f.close()

    sculptor_file=os.path.join(params.directories.temp_dir,"sculpt.ali")
    from phenix.autosol.delete_file import delete_file
    delete_file(sculptor_file)
    f=open(sculptor_file,'w')
    print(align_text, file=f)
    f.close()
    aligned_search_model_folder=params.directories.temp_dir
    aligned_search_model_root="aligned_model"
    aligned_search_model=os.path.join(aligned_search_model_folder,
        aligned_search_model_root+"_sculpt.pdb")
    self.wait_for_file_to_appear(sculptor_file,
         max_ticks=params.control.max_wait_time,allow_failure=False)
    print("Alignment formatted for sculptor: %s\n" %(sculptor_file), file=out)
    args=["model.file_name="+str(pdb_file),
       "alignment.file_name="+str(sculptor_file),
       #"remove_alternate_conformations=True", # add if available 2010-12-03
       "rename=True",
       "renumber.use=original",
       "output.root="+aligned_search_model_root,
       "output.folder="+aligned_search_model_folder,
        ]
    from phaser import sculptor
    from phaser import tbx_utils
    factory = tbx_utils.PhilArgumentFactory( master_phil = sculptor.PHIL_MASTER )
    processed_args=[ factory( argument = arg ) for arg in args ]
    s=sculptor
    file_list=s.run(processed_args,out=out)
    if not file_list or len(file_list)<1 or len(file_list[0])<1:
      raise Sorry("Unable to apply alignment?")
    else:
      print("Aligned model from sculptor: %s" %(str(file_list)), file=out)
    self.wait_for_file_to_appear(aligned_search_model,
         max_ticks=params.control.max_wait_time,allow_failure=False)
    #JUST XXX pass as strings if possible
    new_pdb_input = get_pdb_inp(file_name=aligned_search_model)
    aligned_search_model_text=\
       new_pdb_input.construct_hierarchy().as_pdb_string()
    #JUST
    return aligned_search_model_text

  def decide_to_get_symmetry(self,params,try_to_get_symmetry=None):
    if hasattr(self,'crystal_symmetry') and self.crystal_symmetry:
      return try_to_get_symmetry

    if try_to_get_symmetry is None:
      # if mr_rosetta and not  model_already_placed then do not get symmetry
      if hasattr(params,'place_model') and \
          params.place_model.model_already_placed:
        try_to_get_symmetry=True
    return try_to_get_symmetry

  def start_run(self,run_file):
      from libtbx import easy_run
      easy_run.fully_buffered(
         self.single_run_command+" "+self.add_double_quote(run_file,escape_space=False)).raise_if_errors()

  def as_list(self,value):
    if not value: return []
    r_list=[]
    for item in value:
      r_list.append(item)
    return r_list


  def set_result_file(self,params,out=sys.stdout,stage=None):
    # write to workdir unless it is blank and this is not a sub_process

    if stage is None or params.non_user_params.is_sub_process:
      file='results.pkl'
    else:
      file=stage+'_results.pkl'

    if self.params.directories.output_dir:
      result_file=os.path.join(self.params.directories.workdir,file)
    else:
      result_file=self.make_full_path(file)
    result_csv=result_file[:-4]+".csv"
    return result_file,result_csv

  def write_results(self,params,out=sys.stdout,stage=None):
    # write out solutions as pkl and as csv
    saved_result_file,saved_csv_file=self.set_result_file(params,
       out=out,stage=stage)
    self.save_result(params,saved_result_file=saved_result_file,
      result=self.results)
    self.write_csv(params,saved_csv_file=saved_csv_file,out=out)
    print("\nSaved overall mr_rosetta results in %s" %(saved_result_file), file=out)
    print("\nTo see details of these results type\n"+\
         "    phenix.mr_rosetta mr_rosetta_solutions=%s "%(saved_result_file)+\
       " display_solutions=True\n", file=out)

  def save_result(self,params,saved_result_file=None,result=None):
    from libtbx import easy_pickle
    if params.non_user_params.write_local_files:
      if type(result)==type([1,2,3]):
        for result_indiv in result:
           result_indiv.make_files_local()
      else:
        result.make_files_local()
      easy_pickle.dump(saved_result_file,result)
    else:
      easy_pickle.dump(saved_result_file,result)
    self.result_file=saved_result_file
    libtbx.call_back(message="mr_rosetta_solutions",
      accumulate=False,
      data=os.path.abspath(saved_result_file))

  def setup_work_dir(self,params,out=sys.stdout,prefix="WORK"):
    # Do we need a work directory to work in
    if not params.non_user_params.is_sub_process and  \
          (self.make_full_path(params.directories.workdir)==os.getcwd() or
           not params.directories.workdir):
        params.directories.workdir=self.create_temp_directory(
          n_dir_max=params.control.n_dir_max,
          directory_number=params.control.wizard_directory_number,
          path="",prefix=prefix)
        if not params.directories.top_output_dir or \
         params.directories.top_output_dir==os.getcwd():
          params.directories.top_output_dir=params.directories.workdir

  def setup_temp_dir(self,params,out=sys.stdout,original_directory=None,
     prefix="TEMP"):
    # Do we need a temp directory to work in
    needed=False
    if not params.directories.temp_dir: needed=True
    elif original_directory is not None and  \
           self.make_full_path(params.directories.temp_dir)== \
              self.make_full_path(original_directory): needed=True
    if needed:
      params.directories.temp_dir=self.create_temp_directory(
          n_dir_max=params.control.n_dir_max,
          path="",prefix=prefix)

  def special_cases(self,args,out=sys.stdout):
    # special cases for input files so user doesn't need to specify:
    new_args=[]
    for arg in args:
      if (os.path.isfile(arg)):
        if arg[-3:] in ['sca','mtz']:
          arg_use='input_files.data='+arg
          print("Guessing %s" %(arg_use), file=out)
        elif arg[-3:]=='pdb':
          arg_use='input_files.search_models='+arg
          print("Guessing %s" %(arg_use), file=out)
        elif arg[-3:]=='dat':
          arg_use='input_files.seq_file='+arg
          print("Guessing %s" %(arg_use), file=out)
        else:
          arg_use=arg
      else:
        arg_use=arg
      new_args.append(arg_use)
    return new_args

  def run_subprocesses(self,r=None,out=sys.stdout,title="",
      return_as_groups=False,quiet=False):

    print(78*"=","\nStarting sub-processes %s..." %(title),"\n",78*"=", file=out)
    self.remove_stopwizard()
    try:
      r.run(out=out,quiet=quiet)
    except KeyboardInterrupt:
      self.stop_after_KeyboardInterrupt()
    except Exception as extra_info:
      raise Sorry("mr_rosetta failed...\n"+str(extra_info))

    # now collect all the results
    print(78*"=","\n",\
       "DONE with running subprocesses %s\n" %(title),78*"=", file=out)
    list_of_results=[]
    for result in r.list_of_results:
      if type(result)==type([1,2,3]) and not return_as_groups:
        for rr in result:
         list_of_results.append(rr)
      else:
        list_of_results.append(result)
    return list_of_results

  def print_paths(self,params,out=sys.stdout):
    print("Checking rosetta paths:", file=out)
    print("  rosetta binary: %s " % self.rosetta_path, file=out)
    if hasattr(self,'rosetta_database_dir'):
      print("  database_dir: %s \n  script_dir: %s\n" \
          % ( self.rosetta_database_dir,
             self.rosetta_script_dir), file=out)
    if params.control.rosetta_3_6_or_later:
      print("Running mr_rosetta with inputs for Rosetta 3.6 or later", file=out)

  def find_rosetta(self,params,out=sys.stdout):
    # we need rosetta_path to run...so get it first

    # 2013-10-29 Standard locations of files moved in rosetta 36

    #   rosetta_binary_dir = "rosetta_source/bin" -> main/bin
    #   rosetta_script_dir = rosetta_source/src/apps/public/electron_density -> main/source/src/apps/public/electron_density
    # rosetta_database_dir = "rosetta_database" -> main/database

    if not params.directories.rosetta_path or not \
        os.path.isdir(params.directories.rosetta_path):
      params.directories.rosetta_path=os.getenv("PHENIX_ROSETTA_PATH")
      if not params.directories.rosetta_path or \
         not os.path.isdir(params.directories.rosetta_path):
        raise Sorry("\nSorry, cannot locate rosetta paths.\n"+
         "\nPlease set the environment variable 'PHENIX_ROSETTA_PATH' to indicate "+
         "\nwhere rosetta is to be found. In csh/tcsh use something like:"+
         "\n  setenv PHENIX_ROSETTA_PATH /Users/Shared/unix/rosetta"+
         "\nIn bash/sh use: "+
         "\n  export PHENIX_ROSETTA_PATH=/Users/Shared/unix/rosetta\n")

    if params.control.verbose:
      print("\nRosetta path (params.directories.rosetta_path): %s" %(
         params.directories.rosetta_path), file=out)
    # Now we need to identify the exact rosetta binary. It should be similar
    # to rosetta_binary, but may have different suffixes
    if not params.directories.rosetta_binary_name:
      raise Sorry("Sorry, need to specify the name of rosetta binary ..."+
          "\n(usually 'mr_protocols')")

    target_binary=params.directories.rosetta_binary_name.split(".")[0]+"."
    if params.control.verbose:
       print("Target binary: %s " %(target_binary), file=out)
    rosetta_binary_dir=os.path.split(
       os.path.join(
         params.directories.rosetta_path,params.directories.rosetta_binary_dir,
         params.directories.rosetta_binary_name))[0]

    if not rosetta_binary_dir or not os.path.isdir(rosetta_binary_dir):
      # check out location for version 36 and higher
      old_dir=params.directories.rosetta_binary_dir
      params.directories.rosetta_binary_dir=\
        params.directories.rosetta_binary_dir.replace(
        "rosetta_source",os.path.join("main","source"))
      rosetta_binary_dir_36=os.path.split(
       os.path.join(
         params.directories.rosetta_path,params.directories.rosetta_binary_dir,
         params.directories.rosetta_binary_name))[0]
      if not rosetta_binary_dir_36 or not os.path.isdir(rosetta_binary_dir_36):
        raise Sorry(
         "Sorry...cannot locate the rosetta binary directory as %s or %s"
          %(old_dir,rosetta_binary_dir))
      else:
        rosetta_binary_dir=rosetta_binary_dir_36
        params.control.rosetta_3_6_or_later=True

    if params.control.verbose:
      print("Rosetta binary directory: %s" %(rosetta_binary_dir), file=out)
      if params.control.rosetta_3_6_or_later:
        print("This is Rosetta 3.6 or later", file=out)
      else:
        print("This is pre-Rosetta 3.6", file=out)
    possible_binaries=os.listdir(rosetta_binary_dir)
    rosetta_binary_name=None
    for binary in possible_binaries:
      if binary.find(target_binary)==0:
        rosetta_binary_name=binary
        break
    if not rosetta_binary_name:
      raise Sorry("Sorry...cannot locate a binary starting with '%s' in "
          %(target_binary) +
          "\nthe directory %s " %(rosetta_binary_dir))

    params.directories.rosetta_binary_name=rosetta_binary_name
    if params.control.verbose:
      print("The Rosetta binary name is: %s" %(
        params.directories.rosetta_binary_name), file=out)

    if hasattr(params.directories,'binary_name') and \
        not params.directories.binary_name:
      params.directories.binary_name=rosetta_binary_name.split(".")[-1]
      if params.control.verbose:
        print("Binary name set to: %s" %(params.directories.binary_name), file=out)

    self.rosetta_path=os.path.join(params.directories.rosetta_path,
       params.directories.rosetta_binary_dir,
       params.directories.rosetta_binary_name)

    if params.control.verbose:
      print("Checking on parts of Rosetta path...", file=out)
      base_path=params.directories.rosetta_path
      print("Base path: %s: %s" %(base_path,os.path.exists(base_path)), file=out)
      binary_dir_path=os.path.join(
        base_path,params.directories.rosetta_binary_dir)
      print("Binary directory path: %s: %s" %(
        binary_dir_path,os.path.exists(binary_dir_path)), file=out)
      print("Full Rosetta path: %s: %s\n" %(
        self.rosetta_path,os.path.isfile(self.rosetta_path)), file=out)

    if not os.path.isfile(self.rosetta_path):
      raise Sorry("Sorry...the desired Rosetta binary %s is missing" %(
        self.rosetta_path))

    if hasattr(params.directories,'rosetta_script_dir'):
      self.rosetta_script_dir=os.path.join(params.directories.rosetta_path,
         params.directories.rosetta_script_dir)
      if not os.path.isdir(self.rosetta_script_dir):
        old_dir=params.directories.rosetta_script_dir
        params.directories.rosetta_script_dir=params.directories.rosetta_script_dir.replace( "rosetta_source",os.path.join("main","source"))
        self.rosetta_script_dir=os.path.join(params.directories.rosetta_path,
         params.directories.rosetta_script_dir)
        if not os.path.isdir(self.rosetta_script_dir):
          raise Sorry(
            "Sorry, cannot locate rosetta script directory at %s or %s" %(
             os.path.join(params.directories.rosetta_path,
             params.directories.rosetta_script_dir),
             os.path.join(params.directories.rosetta_path,old_dir)))
        else:
          params.control.rosetta_3_6_or_later=True

    if hasattr(params.directories,'rosetta_database_dir'):
      self.rosetta_database_dir=os.path.join(params.directories.rosetta_path,
         params.directories.rosetta_database_dir)
      if not os.path.isdir(self.rosetta_database_dir):
        old_dir=params.directories.rosetta_database_dir
        params.directories.rosetta_database_dir=params.directories.rosetta_database_dir.replace( "rosetta_database",os.path.join("main","database"))
        self.rosetta_database_dir=os.path.join(params.directories.rosetta_path,
         params.directories.rosetta_database_dir)
        if not os.path.isdir(self.rosetta_database_dir):
          raise Sorry(
         "Sorry, cannot locate rosetta database directory at %s or %s" %(
             os.path.join(params.directories.rosetta_path,
             params.directories.rosetta_database_dir),
             os.path.join(params.directories.rosetta_path,old_dir)))
        else:
          params.control.rosetta_3_6_or_later=True
    # decide if rosetta_3_6_or_later:
    if params.control.rosetta_3_6_or_later is None:
      if params.directories.rosetta_database_dir.find("rosetta_database")>-1:
        params.control.rosetta_3_6_or_later=False
      else:
        params.control.rosetta_3_6_or_later=True
        print("Guessing this is Rosetta 3.6 or later", file=out)

    self.print_paths(params=params,out=out)

  def get_placed_model_list(self,rosetta_solutions):
    placed_model_list=[]
    for solution in rosetta_solutions:
      if not solution.placed_model in placed_model_list:
        placed_model_list.append(solution.placed_model)
    return placed_model_list

  def get_rosetta_solutions_by_group(self,solution_list=[],stage=None,
            score_type=None,group=None,placed_model=None):
    # rosetta solutions by group (chain type)
    solutions_by_group=[]
    for solution in solution_list:
      if solution.group == group and solution.stage==stage \
          and solution.placed_model==placed_model:
        solutions_by_group.append(solution)
    if solutions_by_group:
      sorted_solutions_by_group=self.sort_solutions(
           solution_list=solutions_by_group,
           stage=stage,score_type=score_type)
    else:
      sorted_solutions_by_group=[]
    return sorted_solutions_by_group


  def sort_solutions_by_set(self,solution_list=[],stage="rescored_mr_solution",
      score_type='MR_LLG',any_score_type_if_none=False):
    # rosetta solutions by set (placed_model)

    placed_model_list=self.get_placed_model_list(solution_list)
    solutions_by_set=[]
    for p in placed_model_list:
      solutions_in_set=[]
      for solution in solution_list:
        if solution.placed_model==p:
          solutions_in_set.append(solution)
      sorted_solutions_in_set=self.sort_solutions(
        solution_list=solutions_in_set,stage=stage,
        score_type=score_type,any_score_type_if_none=any_score_type_if_none)
      solutions_by_set.append(sorted_solutions_in_set)
    return solutions_by_set

  def get_best_solution(self,out=None,verbose=False):
    # get best overall Allow score type 2015-10-28

    if self.params.output_files.sort_score_type:
      score_type=self.params.output_files.sort_score_type
    else:
      score_type='R_FACTOR'

    sorted_solution_list=self.sort_solutions(
       stage='ALL',
       solution_list=self.input_solutions+self.results,
       score_type=score_type)
    if sorted_solution_list:
      best_solution=sorted_solution_list[0]
      if best_solution and out is not None:
        print("\nBest solution: ------------------------", file=out)
        self.print_autobuilding_results(result=best_solution,out=out)
        if verbose: best_solution.show_all(out=out)
        print("\n---------------------------------------", file=out)
      return best_solution
    else:
      return None

  def sort_solutions(self,solution_list=[],stage="rescored_mr_solution",
      score_type='MR_LLG',any_score_type_if_none=False):
    # sort solutions based on quality at stage specified
    if not solution_list: return []
    pair_list=[]
    sorted_solutions=[]
    if score_type in ['MR_LLG']:
      reverse=True
    else:
      reverse=False
    any_solutions_at_stage=[]
    for solution in solution_list:
      if stage.lower()=="all" or solution.stage==stage:
        any_solutions_at_stage.append(solution)
        if score_type=='MR_LLG':
           score=solution.score()
        elif score_type=='ROSETTA SCORE':
           score=solution.rosetta_score
        elif score_type == "R_FACTOR" :
          score = solution.r
        else:
           score=solution.score()
        if score is not None:
          pair_list.append([score,solution])
    if pair_list:
      pair_list.sort()
      if reverse:  # sort ascending
        pair_list.reverse()
      for [score,solution] in pair_list:
        sorted_solutions.append(solution)
    elif any_score_type_if_none or stage=="dummy_solution":
       # just get something
       sorted_solutions=any_solutions_at_stage

    return sorted_solutions

  def rewrite_seq_file(self,params,out=sys.stdout):
    # rewrite sequence file because user may have illegal formatting
    seq_file=params.input_files.seq_file
    if not seq_file or not os.path.isfile(seq_file):
      raise Sorry("Sorry the sequence file %s does not seem to exist?"
        %(str(seq_file))  )

    if os.path.split(seq_file)[-1][:7]=="EDITED_":
      return

    if not os.path.isdir(params.directories.temp_dir):
      os.mkdir(params.directories.temp_dir)
    seq_file_name=os.path.split(seq_file)[1]
    new_seq_file=os.path.join(
       params.directories.temp_dir,'EDITED_'+seq_file_name)

    # make it all upper case and fasta format if possible
    # also get rid of any spaces
    # 2011-04-16 also get rid of any "*" as in pir last line
    text=""
    first=True
    for line in open(seq_file).readlines():
      line=line.rstrip()
      if line.replace(" ","")=="": # call it a comment line
        line = "> sequence information"
      elif line[0]==">":  # it is a comment line
        pass
      else:
        line=line.upper().replace(" ","") # 2010-11-18 remove spaces
        line=line.replace("*","") # 2011-04-16 remove "*"
        if first:
          line="> sequence information \n"+line  # add header line
      first=False
      text+=line+"\n"
    f=open(new_seq_file,'w')
    print(text, file=f)
    f.close()
    print("Sequence rewritten to %s : \n %s "%(new_seq_file,text), file=out)
    params.input_files.seq_file=new_seq_file

  def get_file_and_score_list(self,params,temp_dir="",
     nstruct=1,out=sys.stdout):
    id_list=None
    file_and_score_list=[]
    scoring_file=os.path.join(temp_dir,'score.sc')
    self.wait_for_file_to_appear(scoring_file,
       max_ticks=params.control.max_wait_time,allow_failure=True)

    if not os.path.isfile(scoring_file):
      return file_and_score_list  # give up

    for line in open(scoring_file).readlines():
      if not line or len(line.split())<2: continue
      if line.split()[0]=="SEQUENCE:": continue
      spl=line.split()
      if spl[0]!="SCORE:":
        raise Sorry("Unexpected contents of scoring file:\nFILE: %s \n" %(
          scoring_file) + "LINE: %s" %(line) )
      if not id_list:
        id_list=spl[1:]
      else:
        score=None
        description=None
        for id,value in zip(id_list,spl[1:]):
          if id=="total_score":
            score=float(value)
          if id=="description":
            description=value
        if id is not None and description is not None:
          file=os.path.join(temp_dir,description+".pdb")
          if params.non_user_params.dummy_rosetta:
            dummy_file=params.input_files.model
            print("Using %s instead of %s for scoring" %(dummy_file,file), file=out)
            file=dummy_file
          file_and_score_list.append([file,score])

    return file_and_score_list

  def check_fragment_files(self,params,run_rebuild=True,out=sys.stdout):
      if len( params.input_files.fragment_files)<2:
        if params.input_files.fragment_files_chain_list:
           if not params.input_files.fragment_files_3_mer_by_chain or \
               not params.input_files.fragment_files_9_mer_by_chain:
             raise Sorry("Sorry, if you specify "+
             "params.input_files.fragment_files_chain_list you also "+
             "need fragment_files_3_mer_by_chain and "+
             "fragment_files_9_mer_by_chain")
           elif len(params.input_files.fragment_files_chain_list) != \
                len(params.input_files.fragment_files_3_mer_by_chain) or \
                len(params.input_files.fragment_files_chain_list) != \
                len(params.input_files.fragment_files_9_mer_by_chain):
             raise Sorry("Sorry, if you specify "+
             "params.input_files.fragment_files_chain_list you need the same "+
             "numbers of entries in fragment_files_3_mer_by_chain and "+
             "fragment_files_9_mer_by_chain")
           else:
             print("Using fragments specific for each chain")
             return

        if params.control.rosetta_3_6_or_later:
          print("\nNote: fragment files will be generated by Rosetta\n",\
             "(Requires version 2013wk35 or later of Rosetta)\n", file=out)
          params.control.generate_fragment_files=True
          return
        elif params.control.generate_fragment_files:
          raise Sorry("Sorry, you need rosetta 3.6 or later (2013wk35) for "+
            " generate_fragment_files")

        elif (not params.input_files.use_dummy_fragment_files) and \
          hasattr(self,'alignment_file_needed') and self.alignment_file_needed:
          if not params.rosetta_rebuild.run_rosetta_rebuild  or \
             not run_rebuild:
            return
          if hasattr(self,'will_be_run') and \
              not self.will_be_run(params,test_point='run_rebuild'):
             return
          raise Sorry("Sorry, need 2 fragment files (9 and 3 residues) "+
             "for rosetta_rebuild if sequence is not identical to template")
        else:  # just use dummy ones
          import libtbx
          params.input_files.fragment_files=[]
          for k in [3,9]:
            file=libtbx.env.find_in_repositories(
             relative_path=
               "phenix/phenix/rosetta/dummy_fragments_file_%d.gz" %(k),
               test=os.path.isfile)
            params.input_files.fragment_files.append(file)
          print("\nNote: dummy fragments files will be used: ", \
            params.input_files.fragment_files, file=out)
          print("(This is ok if your template has no gaps)", file=out)


      params.input_files.fragment_files=\
      self.make_full_paths_in_list_and_check( self.as_list(
          params.input_files.fragment_files),'fragment_files')
      # now make sure they are in order...big first then small
      # NOTE: must have run make_full_paths to make this a list first
      if params.input_files.sort_fragment_files:
        params.input_files.fragment_files.sort()
        params.input_files.fragment_files.reverse()


  def simple_get_model_info(self,params,out=sys.stdout,
      try_to_get_symmetry=None):

      # create model_info_object from our PDB, alignment and seq file
      print("\nCreating information file ", file=out)
      hierarchy = self.pdb_input.construct_hierarchy()
      crystal_symmetry=self.none_if_no_space_group_info(
        self.pdb_input.crystal_symmetry_from_cryst1())
      if crystal_symmetry:
             if try_to_get_symmetry:
               print("SETTING CRYSTAL SYMMETRY FROM MODEL",\
                iotbx.pdb.format_cryst1_record(
                crystal_symmetry=crystal_symmetry), file=out)
               self.crystal_symmetry=crystal_symmetry

      self.model_info=self.create_dict_from_hierarchy(hierarchy=hierarchy,out=out)
      if not hasattr(self,'crystal_symmetry'):
        self.crystal_symmetry=None

      self.set_model_info()

  def edit_rosetta_model_from_solution(self,params,solution=None,
      out=sys.stdout,add_id=None):
    if os.path.split(solution.model)[-1][-7:] == "_ed.pdb":
      return solution.model # already edited

    if add_id:
      id=solution.id
    else:
      id=None
    rewritten_model=self.apply_model_info(params,model=solution.model,
      apply_ncs=True,ncs_object=solution.ncs_object,out=out,
      crystal_symmetry=solution.crystal_symmetry,
      set_chain_as_A=True,id=id)  # 2010-12-16 set one_chain_only
    return rewritten_model

  def get_crystal_symmetry_from_pdb(self,model=None,pdb_input=None,
      out=sys.stdout,require_symmetry=True):
      self.crystal_symmetry=self.none_if_no_space_group_info(
        pdb_input.crystal_symmetry_from_cryst1())
      if self.crystal_symmetry is None:
        if require_symmetry:
           raise Sorry("Need crystal symmetry in"+
          " %s for get_crystal_symmetry_from_pdb" %(str(model)))
      else:
        print("\nSETTING CRYSTAL SYMMETRY FROM %s" %(model), file=out)
        print(iotbx.pdb.format_cryst1_record(
           crystal_symmetry=self.crystal_symmetry), file=out)

  def add_ed(self,model,ext="_ed.pdb",id=None):
       from phenix.autosol.trim_file_name import trim_file_name
       trimmed=trim_file_name(model).trimmed_file
       if id is not None:
         ext="_%s%s" %(str(id),ext) # add id
       if len(trimmed)>8:  # need to keep at least 5 chars identical
         rewritten_model=model[:-4]+ext
       else:
         rewritten_model=model+ext
       return rewritten_model

  def apply_model_info(self,params,model=None,rewritten_model=None,
      apply_ncs=False,ncs_object=None,set_chain_as_A=False,out=sys.stdout,
      crystal_symmetry=None,id=None):
    if rewritten_model is None:
       rewritten_model=self.add_ed(model,id=id)
    # we are going to copy something (in this case atom.b) from model_info to
    # model to yield rewritten_model
    # if no model_info just copy
    # get rid of hydrogens if possible

    # As current model may have different numbering from original model, use
    # the dictionary orig_resno_from_aligned_resno_dict to get original resno

    print("\nEditing %s and writing to %s\n" %(
       model,rewritten_model), file=out)
    pdb_input = get_pdb_inp(file_name=model)
    hierarchy = pdb_input.construct_hierarchy()

    # chain:id
    # residue_group: resseq icode
    # atom_group: altloc resname
    # atom: xyz occ b segid format_atom_record()
    n_mean=0
    n_used=0
    f=open(rewritten_model,'w')
    if not crystal_symmetry:
      crystal_symmetry=self.dummy_crystal_symmetry()
    print(iotbx.pdb.format_cryst1_record(
         crystal_symmetry=crystal_symmetry), file=f)
    print("CRYSTAL SYMMETRY: ",iotbx.pdb.format_cryst1_record(
         crystal_symmetry=crystal_symmetry))

    for model in hierarchy.models()[:1]:
        for chain in model.chains():
          if set_chain_as_A or not chain.id:
            chain.id="A"  # setting chain ID
          for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
              for atom in atom_group.atoms():
                if atom.element_is_hydrogen():  # doesn't work 2010-08-04 # now it does 2013-11-15
                   continue
                i=residue_group.resseq_as_int()
                original_resseq_int=\
                  self.orig_resno_from_aligned_resno_dict.get(i,-999)
                id=self.get_atom_id(
                  chain_id=chain.id,
                  resseq_int=original_resseq_int,
                  #icode=residue_group.icode,
                  #resname=atom_group.resname,
                  #altloc=atom_group.altloc,
                  name=atom.name)
                if id in self.model_info.keys():
                  atom.b=self.model_info[id]
                  n_used+=1
                elif 'mean' in self.model_info.keys():
                  atom.b=self.model_info['mean']
                  n_mean+=1
                else:
                  raise Sorry("Sorry, model info dict must have 'mean' as a key")
                if not self.looks_like_hydrogen(atom.name):
                  print(atom.format_atom_record(), file=f)
    f.close()
    self.wait_for_file_to_appear(rewritten_model,
       max_ticks=params.control.max_wait_time,allow_failure=False)
    print("Total of %d atoms matched to info_file and %d not" %(
         n_used,n_mean), file=out)

    # NOW APPLY NCS IF PRESENT TO GENERATE ENTIRE MOLECULE 2010-09-06
    if apply_ncs and ncs_object: # overwrite rewritten model
      self.apply_ncs_to_molecule(params,rewritten_model=rewritten_model,
        ncs_object=ncs_object,out=out)

    # ZZZ 2013-09-12 REMOVE BAD CONTACTS IF REQUESTED
    return rewritten_model


  def apply_ncs_to_molecule(self,params,rewritten_model=None,
        ncs_object=None,out=sys.stdout):
    pdb_input = get_pdb_inp(file_name=rewritten_model)
    from phenix.autosol.delete_file import delete_file
    delete_file(rewritten_model)
    copies_placed=self.generate_ncs_copies(params,pdb_input=pdb_input,
      output_pdb=rewritten_model,
      ncs_object=ncs_object,out=out)
    print("\nFile with all NCS copies in %s " %(rewritten_model), file=out)
    return rewritten_model

  def generate_ncs_copies(self,params,pdb_input=None,output_pdb=None,
    ncs_object=None,out=sys.stdout,max_copies=None,start_copy=None):
    # create all the copies from one copy
    from phenix.command_line.apply_ncs import apply_ncs
    args=["temp_dir="+params.directories.temp_dir,
          "pdb_out="+output_pdb]
    if max_copies is not None:
       args.append("max_copies="+str(max_copies))
    if start_copy is not None:
       args.append("start_copy="+str(start_copy))
    an=apply_ncs(args,ncs_object=ncs_object,pdb_input=pdb_input,out=out)
    # XXX NOTE: need to apply_ncs using correct ncs group if multiple chain types 2010-12-11
    self.wait_for_file_to_appear(output_pdb,
         max_ticks=params.control.max_wait_time,allow_failure=False)
    return an.copies_placed

  def looks_like_hydrogen(self,text): # H Hanything 123456789H are all H
    t=text.rstrip().upper().replace(" ","")
    if t[0]=='H':
      return True
    elif len(t)>=2:
      if t[1]=='H' and t[0] in "0123456789":
         return True
    return False

  def set_model_info(self):
    from phenix.rosetta.mr_rosetta import model_info_object
    self.model_info_object=model_info_object()
    for x in self.model_info_object.headers:
      if hasattr(self,x):
       setattr(self.model_info_object,x,getattr(self,x))

  def create_dict_from_hierarchy(self,hierarchy=None,out=sys.stdout):
    # create a dictionary of hashed atom identifiers so we can find info on
    # any atom quickly. Save only B factors for now...but can save anything

    # to use....get ID for an atom by calling get_atom_id and then you have result.
    # if no result, use 'mean'....

    dd={}
    residue_list=[]
    mean=0.
    mean_n=0.
    for model in hierarchy.models()[:1]:
      for chain in model.chains():
        chain_id=chain.id
        if not chain_id: chain_id="A"
        for residue_group in chain.residue_groups():
          for atom_group in residue_group.atom_groups():
            for atom in atom_group.atoms():
              id=self.get_atom_id(
                 chain_id=chain_id,
                 resseq_int=residue_group.resseq_as_int(),
                 #icode=residue_group.icode,
                 #resname=atom_group.resname,
                 #altloc=atom_group.altloc,
                 name=atom.name)
              #print >>out,"wID:",id,atom.format_atom_record()
              info=atom.b # put whatever you want here
              mean+=atom.b
              mean_n+=1.
              dd[id]=info
    if mean_n>1: mean=mean/mean_n
    dd['mean']=mean
    print("\nMean B-value for this structure was %6.1f " %(mean), file=out)
    return dd

  def get_atom_id(self,chain_id="",
       resseq_int=0,
       icode="",
       resname="",
       altloc="",
       name=""):
    # residue_group: resseq icode
    # atom_group: altloc resname
    # atom: xyz occ b segid format_atom_record()
    if chain_id.replace(" ","")=="": chain_id="A"
    atom_id=":".join([str(chain_id).replace(" ",""),
                      str(resseq_int).replace(" ",""),
                      str(icode).replace(" ",""),
                      str(resname).replace(" ",""),
                      str(altloc).replace(" ",""),
                      str(name).replace(" ","")])
    return atom_id

  def get_res_id(self,chain_id,resseq):
     return ":".join([str(chain_id),str(resseq)])


  def create_dummy_alignment_lines(self,base_name='coords',
      target_sequence=None,pdb_sequence=None):
    all_lines=[]
    if target_sequence and pdb_sequence:
      target_sequence_list=[target_sequence]
      pdb_sequence_list=[pdb_sequence]
    else:
      target_sequence_list=self.sequence
      pdb_sequence_list=self.sequence

    for target_sequence,pdb_sequence in zip(
          target_sequence_list,pdb_sequence_list):
        # 2011-01-06 check to make sure no gaps
        all_lines.append("## TARGET %s \n" %(base_name))
        all_lines.append("# hhsearch\n")
        all_lines.append("scores_from_program: 0 1.00\n")
        all_lines.append("0 %s\n" %(target_sequence))
        all_lines.append("0 %s\n" %(pdb_sequence))
        all_lines.append("--\n")
    return all_lines

  def read_seq_file(self,params,out=sys.stdout,quiet=True):
    if not params.input_files.seq_file or \
       not os.path.isfile(params.input_files.seq_file):
      raise Sorry("\nSorry, the sequence file '%s' is missing?" %(
        params.input_files.seq_file))

    if not quiet:
      print(80*"=", file=out)
      print("\nReading sequence from %s" %(params.input_files.seq_file), file=out)
      print(80*"=", file=out)
      print("\nSequence info: \n%s" %(
         open(params.input_files.seq_file).read()), file=out)
    self.Facts={'verbose':params.control.verbose,
                'chain_type':'PROTEIN'}
    chains=self.Build_read_seq(seq_file=params.input_files.seq_file,
      set_start_chains=False,out=out)
    if not chains:
      raise Sorry("\nSorry, there does not appear to be any sequence in %s?" %(
       params.input_files.seq_file))
    #print "SEQ: ",chains
    if len(chains)>1:
      print("NOTE: Multiple chains in "+\
        " %s: %s" %(params.input_files.seq_file,chains), file=out)
    self.sequence=chains


  def simple_read_alignment_file_and_model(self,params,alignment_file=None,
     model=None,out=sys.stdout,quiet=True, try_to_get_symmetry=None):

    self.simple_check_alignment_file_read_model(params,
         alignment_file=alignment_file,model=model,quiet=quiet,out=out)
    # now we have set self.sequence_by_group etc...
    # also self.pdb_input
    identity=self.check_seq_vs_pdb(out=out,quiet=quiet) # checks seq vs PDB
    if hasattr(params,'place_model') and params.place_model.identity is None:
      params.place_model.identity=identity
      print("Identity is: %d percent" %(identity), file=out)

    self.simple_get_model_info(params,out=out,
        try_to_get_symmetry=try_to_get_symmetry)

  def check_seq_vs_pdb(self,out=sys.stdout,quiet=True):
    # now get new residue numbers (and residue types) for each residue in
    # our PDB file (after alignment)
    pdb_chain_dict_keys=self.pdb_chain_dict.keys()
    pdb_chain_dict_keys.sort()
    self.orig_resno_from_aligned_resno_dict={}
    pdb_residue_list_by_chain_id=\
       self.check_pdb_chain_dict_and_set_chain_id()
    self.get_group_by_chain_id()
    for chain_id in self.pdb_chain_dict.keys():
      # figure out which group this is associated with
      group=self.group_by_chain_id[chain_id]
      sequence=self.sequence_by_group[group]
      seq_file_start=self.seq_file_start_by_group[group]
      pdb_file_start=self.pdb_file_start_by_group[group]
      pdb_file_alignment=self.pdb_file_alignment_by_group[group]
      seq_file_alignment=self.seq_file_alignment_by_group[group]

      pdb_residue_list=self.pdb_chain_dict[chain_id][1]
      pdb_resno_list=self.pdb_chain_dict[chain_id][0]

      if not quiet:
        print("Getting new residue names and numbers for chain %s" %(chain_id), file=out)
        print("Alignment starts at position %d in sequence file (%s)" %(
        seq_file_start+1,sequence[seq_file_start]), file=out)
        print("...and at position %d in PDB file (residue %s %s)" %(
         pdb_file_start+1,pdb_residue_list[pdb_file_start],
            pdb_resno_list[pdb_file_start]), file=out)

      # Make sure PDB sequences match:
      seq_to_match_from_pdb="".join(pdb_residue_list[pdb_file_start:]).replace("-","")
      seq_to_match_from_align_file=pdb_file_alignment.replace("-","")
      if seq_to_match_from_pdb != seq_to_match_from_align_file:
        print("Trying to adjust the PDB sequence from alignment file:"+\
           "\n%s\nto match the actual PDB sequence:\n%s" %(
            seq_to_match_from_align_file,seq_to_match_from_pdb), file=out)

        #new alignment line will have same number of chars as previous.
        new_alignment_line,matches=self.adjust_alignment_to_match_pdb(
              pdb_file_alignment,seq_to_match_from_pdb,out=out)
        if new_alignment_line:
          pdb_file_alignment=new_alignment_line
        else:
          raise Sorry("\nSorry, the sequence from the PDB file starting at "+
            "position %d \n(%s) " %(pdb_file_start+1,
            seq_to_match_from_pdb)+"\ndoes not match the sequence from "+
            "the alignment file line 5 \n(%s) " %(seq_to_match_from_align_file))

      i=seq_file_start-1
      count_tot=0
      count_same=0
      count_diff=0
      for pdb_file_resno,pdb_file_residue,pdb_file_alignment_residue,\
         seq_file_alignment_residue in zip(
             pdb_resno_list[pdb_file_start:],
             pdb_residue_list[pdb_file_start:],
             pdb_file_alignment,
             seq_file_alignment):
        count_tot+=1
        if pdb_file_alignment_residue=="-":
           # nothing for this residue in PDB file
           i+=1  # residue starts at 1 if seq_file_start=0
        elif seq_file_alignment_residue=="-": # removing this pdb residue
           pass
        else: # both present:
           i+=1  # residue starts at 1 if seq_file_start=0
           self.orig_resno_from_aligned_resno_dict[i]=pdb_file_resno #XXXXHERE
           if pdb_file_alignment_residue==seq_file_alignment_residue:
             count_same+=1
           else:
             count_diff+=1
    if count_diff+count_same:
      identity=int(0.5+100.*float(count_same)/float(count_diff+count_same))
    else:
      identity=0
    if count_tot != count_same:
      self.alignment_file_needed=True

    # list results
    for chain_id in self.pdb_chain_dict.keys():
      # figure out which group this is associated with
      group=self.group_by_chain_id[chain_id]
      sequence=self.sequence_by_group[group]
      seq_file_start=self.seq_file_start_by_group[group]
      pdb_file_start=self.pdb_file_start_by_group[group]
      pdb_file_alignment=self.pdb_file_alignment_by_group[group]
      seq_file_alignment=self.seq_file_alignment_by_group[group]
      pdb_residue_list=self.pdb_chain_dict[chain_id][1]
      pdb_resno_list=self.pdb_chain_dict[chain_id][0]

      for i in range(seq_file_start,
         seq_file_start+len(seq_file_alignment)-seq_file_alignment.count("-")):
        orig_resno=self.orig_resno_from_aligned_resno_dict.get(i,-999)
        orig_res_id=self.get_res_id(chain_id,orig_resno)
        orig_residue=self.res_dict.get(orig_res_id)
        if i+1> len(sequence):
          raise Sorry("\nSorry, sequence, alignment file, and offsets do "+
          "not match\nAlignment sequence for "+
          "sequence_file (line 4 of alignment file) "+
          "\nhas %d residues and %d missing residues" %(
           len(seq_file_alignment)-seq_file_alignment.count("-"),
                 seq_file_alignment.count("-")) +
          "\n(%s)"%(seq_file_alignment) +
          "\nwhile the sequence file itself has only %d residues " %(
             len(sequence)-seq_file_start) +
            "\n(excluding residues before %d) " %(seq_file_start+1) +
           "\n(%s)" %(sequence[seq_file_start:]))
        #print "RESIDUE %s (%s): In PDB: %s (%s)" % (
          #str(i+1),sequence[i],str(orig_resno),str(orig_residue))

    return identity

  def read_model(self,model,params=None,out=sys.stdout,quiet=True,
     one_chain_only=False,try_to_get_symmetry=None):
    # read in model and get sequence vs residue numbers in each chain
    print(80*"=", file=out)
    print("Reading model from %s " %(model), file=out)
    print(80*"=", file=out)
    if not os.path.isfile(model):
      raise Sorry("Sorry, the model %s is missing?" %(model))
    lines=open(model).readlines()
    if hasattr(params,'place_model') and \
         hasattr(params.place_model,'fixed_model') and \
         params.place_model.fixed_model:
       if lines[-1].rstrip()=='END': lines=lines[:-1]
       new_lines=open(params.place_model.fixed_model).readlines()
       for line in new_lines:
        if not line.startswith("CRYST1"): lines.append(line)
    pdb_input = get_pdb_inp(lines=lines)
    self.phaser_identity=self.get_phaser_identity_remark(pdb_input)
    if self.phaser_identity is not None:
       print("PHASER IDENTITY: ",self.phaser_identity)
    hierarchy = pdb_input.construct_hierarchy()
    self.pdb_input=pdb_input
    if try_to_get_symmetry:
      self.get_crystal_symmetry_from_pdb(model=model,pdb_input=pdb_input,
         out=out,require_symmetry=True)

    self.res_dict={}
    # list of chains : [pdb_resseq_list,pdb_resname_list]
    #  2013-02-14 if fixed_model...use that too

    self.pdb_chain_dict={}
    from phenix.autosol.resname import resname
    for model in hierarchy.models()[:1]:
      for chain in model.chains():
        chain_id=chain.id
        if chain_id.replace(" ","")=="": chain_id="A"
        if chain_id in self.pdb_chain_dict.keys(): # add on
          pdb_resseq_list=self.pdb_chain_dict[chain_id][0]
          pdb_resname_list=self.pdb_chain_dict[chain_id][1]
        else:
          pdb_resseq_list=[]
          pdb_resname_list=[]
        self.pdb_chain_dict[chain_id]=[pdb_resseq_list,pdb_resname_list]
        skipped_residues=[]
        for conformer in chain.conformers()[:1]:
          for residue in conformer.residues():
            resname_one=resname(chain_type='PROTEIN',
              semet=False).OneCharName(residue.resname)
            if not resname_one:
              skipped_residues.append("%s %s %s"
                  %(chain_id,residue.resseq,residue.resname))
              continue
            resseq=residue.resseq_as_int()
            pdb_resseq_list.append(resseq)
            pdb_resname_list.append(resname_one)
            res_id=self.get_res_id(chain_id,resseq)
            self.res_dict[res_id]=resname_one
        if skipped_residues and params is not None and params.control.verbose:
          skip_file=os.path.join(
              params.directories.temp_dir,'SKIPPED_RESIDUES.dat')
          ff=open(skip_file,'w')
          print("Total of %d skipped residues and waters listed in: %s \n"%(
                 len(skipped_residues),skip_file), file=out)
          print("Skipped residues", file=ff)
          print(skipped_residues, file=ff)
          ff.close()

    for key in self.pdb_chain_dict.keys():
      if self.pdb_chain_dict[key]==[[],[]]:
        del self.pdb_chain_dict[key]

  def get_phaser_identity_remark(self,pdb_input):
    for line in pdb_input.remark_section():
      if not line.find(' PHASER ENSEMBLE ')>-1:continue
      if not line.find(' ID ')>-1:continue
      spl=line.split(" ID ")
      if len(spl)<2: continue
      id_line=spl[1]
      id=id_line.split()[0]
      try:
       identity=float(id)
      except Exception: continue
      return identity
    return None

  def get_file_base(self,file_name):
      file_split=os.path.split(file_name)[-1].split(".")
      if len(file_split)>1:
        file_base=".".join(file_split[:-1])
      else:
        file_base=file_split[0]
      return file_base

  def get_group_by_chain_id(self,out=sys.stdout):
    self.group_by_chain_id={}
    for chain_id in self.pdb_chain_dict.keys():
      self.group_by_chain_id[chain_id]=self.get_group_from_chain_id(chain_id,
       out=out)


  def make_dict(self,pair):
    # make a dictionary out of a matched pair of lists
    assert len(pair) == 2
    dd={}
    for key,value in zip(pair[0],pair[1]):
      dd[key]=value
    return dd

  def check_pdb_chain_dict_and_set_chain_id(self):
         # require that all the chains have the same residues (no mismatches)
         # as a standard chain.  Return list of unique chain_id values

         pdb_chain_dict_keys=self.pdb_chain_dict.keys()
         pdb_chain_dict_keys.sort()
         chain_id_by_chain_id=[pdb_chain_dict_keys[0]]
         if len(pdb_chain_dict_keys)!=1:
           standard_dict_list=[self.make_dict(
             self.pdb_chain_dict[pdb_chain_dict_keys[0]])]

           for key in pdb_chain_dict_keys[1:]:
             dd=self.make_dict(self.pdb_chain_dict[key])
             is_old=False
             for standard_dict in standard_dict_list:
               in_this_standard_dict=True
               for id in dd.keys():
                 if not id in standard_dict.keys(): continue
                 if standard_dict[id] != dd[id]:
                    in_this_standard_dict=False
               if in_this_standard_dict:
                 is_old=True # found already
             if not is_old:  # this is a new one
               standard_dict_list.append(
                 self.make_dict(self.pdb_chain_dict[key]))
               chain_id_by_chain_id.append(key)


         pdb_residue_list_by_chain_id={}
         for chain_id in chain_id_by_chain_id:
           pdb_residue_list_by_chain_id[chain_id]=\
              self.pdb_chain_dict[chain_id][1]
         return pdb_residue_list_by_chain_id

  def simple_check_alignment_file_read_model(self,params,model=None,
     alignment_file=None,quiet=True,out=sys.stdout):

    # 2011-05-14 allow simple alignment file as used by sculptor...

    # 2010-12-09 allow several alignments in one file
    if alignment_file is not None and not os.path.isfile(alignment_file):
      raise Sorry("Sorry, the alignment file %s is missing? " %(alignment_file))


    if not quiet:
       print(80*"=", file=out)
       print("\nReading alignment from %s" %(str(alignment_file)), file=out)
       print(80*"=", file=out)


    pdb_residue_list=None
    if hasattr(params,'place_model') and \
          params.place_model.model_already_placed:
        try_to_get_symmetry=True
    else:
        try_to_get_symmetry=False

    self.read_model(model,params=params,out=out,quiet=True,
            try_to_get_symmetry=try_to_get_symmetry)

    if params.place_model.identity is None and self.phaser_identity is not None:
      print("Setting identity to %s from phaser REMARK record\n" %(
        str(self.phaser_identity)))
      params.place_model.identity=self.phaser_identity

    pdb_residue_list_by_chain_id=\
            self.check_pdb_chain_dict_and_set_chain_id()
    # read it and edit it in a simple way if we can

    if alignment_file is not None:
      all_lines=open(alignment_file).readlines()
      if len(all_lines)<6 or all_lines[0][:3] != "## ":
        all_lines=self.get_alignment_lines_from_sculptor_file(all_lines,model)
        if not all_lines:
          raise Sorry("\nSorry, the alignment file %s \n" %(
           params.input_files.alignment_files) +
           "does not seem to have the correct format?\n"+
           "It can either start with '##' or with '> title of target seq here'")

    else:  # create dummy lines and dummy file
      all_lines=self.create_dummy_alignment_lines(
         base_name=self.get_solution_key(model))
      dummy_align_file=os.path.join(params.directories.temp_dir,'alignment.ali')
      f=open(dummy_align_file,'w')
      print("".join(all_lines), file=f)
      f.close()
      self.wait_for_file_to_appear(dummy_align_file,
          max_ticks=params.control.max_wait_time,allow_failure=False)
      print("Wrote dummy alignment to %s" %(dummy_align_file), file=out)
      params.input_files.alignment_files=[dummy_align_file]

    groups=len(all_lines)//6
    if groups < 1:
      raise Sorry("\nSorry, the alignment file %s \n" %(
          params.input_files.alignment_files) +
       "does not have enough lines (6 required, %d found):" %(len(all_lines))+
       "\n"+''.join(all_lines)+"\n\nNote: last line needs to be '--'\n")
    matched_sequence_by_group=[]
    self.seq_file_by_group=[]
    self.align_file_by_group=[]

    self.seq_file_start_by_group=[]
    self.pdb_file_start_by_group=[]
    self.seq_file_alignment_by_group=[]
    self.pdb_file_alignment_by_group=[]

    # now save alignment info for all groups
    if params.control.verbose:
      print("\nFinal alignment info:\n%s " %("\n".join(all_lines)), file=out)

    all_text=""
    for group in range(groups):
      alignment_lines=all_lines[group*6:group*6+6]
      new_alignment_lines,matched_sequence=\
          self.simple_interpret_alignment_lines(params,
            alignment_lines=alignment_lines,
            alignment_file=alignment_file,out=out)
      matched_sequence_by_group.append(matched_sequence)
      self.seq_file_start_by_group.append(int(alignment_lines[3].split()[0]))
      self.pdb_file_start_by_group.append(int(alignment_lines[4].split()[0]))
      self.seq_file_alignment_by_group.append(alignment_lines[3].split()[1])
      self.pdb_file_alignment_by_group.append(alignment_lines[4].split()[1])

      if 1: # XXX 2011-01-06
        base_file=os.path.join(params.directories.temp_dir,
         'group_'+str(group))
        seq_file=base_file+".seq"
        align_file=base_file+".ali"

        ff=open(seq_file,'w')
        print(">group "+str(group), file=ff)
        print(matched_sequence, file=ff)
        ff.close()

        ff=open(align_file,'w')
        for line in new_alignment_lines:
          print(line.rstrip(), file=ff)
          all_text+=line.rstrip()+"\n"
        ff.close()
        self.seq_file_by_group.append(seq_file)
        self.align_file_by_group.append(align_file)
        self.wait_for_file_to_appear(align_file,
          max_ticks=params.control.max_wait_time,allow_failure=False)

    # Now write out the final edited alignment file
    edited_align_file=os.path.join(
          params.directories.temp_dir,'edited_align.ali')
    f=open(edited_align_file,'w')
    print("".join(all_text), file=f)
    f.close()
    self.wait_for_file_to_appear(edited_align_file,
          max_ticks=params.control.max_wait_time,allow_failure=False)
    print("\nFully edited alignment file: %s " %(edited_align_file), file=out)
    params.input_files.alignment_files=[edited_align_file]
    print("Matched sequences: ",matched_sequence_by_group, file=out)
    self.sequence_by_group=matched_sequence_by_group
    self.number_of_groups=len(self.sequence_by_group)
    # Now get best match of group for each chain_id
    self.get_group_by_chain_id(out=out)


  def get_alignment_lines_from_sculptor_file(self,all_lines,model):
    # try to interpret lines
    text=""
    for line in all_lines:
      if 0 and line[:1] == "-": # 2011-08-26 this causes problems if line starts with "----" so skip it
        continue
      elif line[:1] == '>':
        text+=">"
      else:
        text+=line.rstrip()
    lines=text.split(">")
    if len(lines) != 3: return ""
    all_lines=self.create_dummy_alignment_lines(
         base_name=self.get_solution_key(model),
         target_sequence=lines[1],pdb_sequence=lines[2])
    return all_lines

  def get_group_from_chain_id(self,chain_id,out=sys.stdout):
    # requires read_alignment_file first
    group=-1
    pdb_residue_list=self.pdb_chain_dict[chain_id][1]
    best_group=None
    best_matches=0
    for pdb_file_alignment,pdb_file_start in zip(
        self.pdb_file_alignment_by_group,self.pdb_file_start_by_group):
      group+=1
      seq_to_match_from_pdb=\
        "".join(pdb_residue_list[pdb_file_start:]).replace("-","")
      new_alignment_line,matches=self.adjust_alignment_to_match_pdb(
              pdb_file_alignment,seq_to_match_from_pdb,out=out,printout=False)
      if matches>best_matches:
        best_matches=matches
        best_group=group
    return best_group

  def get_best_pdb_residue_list(self,pdb_file_alignment,pdb_file_start,
     out=sys.stdout):
    # return best-matching chain from pdb file
    best_seq_to_match_from_pdb=None
    best_matches=0
    best_new_alignment_line=None
    seq_to_match_from_pdb="" # 2012-04-03 in case we have no chains..
    for chain_id in self.pdb_chain_dict.keys():
      pdb_residue_list=self.pdb_chain_dict[chain_id][1]
      seq_to_match_from_pdb=\
        "".join(pdb_residue_list[pdb_file_start:]).replace("-","")


      new_alignment_line,matches=self.adjust_alignment_to_match_pdb(
              pdb_file_alignment,seq_to_match_from_pdb,out=out,printout=False)
      if seq_to_match_from_pdb.replace("-","")==\
           pdb_file_alignment.replace("-",""):
         print("Alignment file matches PDB file", file=out)
         return seq_to_match_from_pdb,None
      elif matches > best_matches:
         best_new_alignment_line=new_alignment_line
         best_matches=matches
         best_seq_to_match_from_pdb=seq_to_match_from_pdb

    print("Trying to adjust the PDB sequence from alignment file:"+\
           "\n%s\nto match the actual PDB sequence:\n%s" %(
            pdb_file_alignment,seq_to_match_from_pdb), file=out)
    if not best_new_alignment_line:
      raise Sorry("Sorry line 5 of the alignment file has the sequence \n%s" %(
        pdb_file_alignment) + "\nwhich does not match the sequence from the "+\
        "pdb file: \n%s" %(seq_to_match_from_pdb))
    return best_seq_to_match_from_pdb,best_new_alignment_line

  def simple_interpret_alignment_lines(self,params,alignment_lines=None,
        alignment_file=None,out=sys.stdout):

    # try to read alignment file just as rosetta does. This gives us
    # the residue number correspondences for apply_model_info...
    # Looks like:
    standard_align="""
## RB_ coords.pdb
# hhsearch
scores_from_program: 0 1.00
1 QLMDMRD
0 QLMDMRD
--
    """
    # and goes with sequence file like:
    standard_seq="""
> coords fasta
AQLMDMRD
    """

    aa=[]
    for line in alignment_lines:
      if line.rstrip(): aa.append(line)
    alignment_lines=aa
    print("ALIGNMENT LINES: \n%s" %("".join(aa)))

    # identify which sequence this matches
    spl,start_pos,matched_sequence=self.get_matched_sequence(alignment_lines)
    if start_pos<0:
        raise Sorry("Sorry, the alignment file %s has a target sequence \n(%s)"\
          %(alignment_file,alignment_lines[3]) +"\n that does not match the "+\
          "sequence in the sequence file \n(%s)" %(str(self.sequence)))
    if len(alignment_lines[3].split()) != 2 or len(alignment_lines[4].split()) != 2:
      raise Sorry("\nSorry, the seq_file and pdb_file alignment lines \n"+
       " %s\n%s " %(alignment_lines[3],alignment_lines[4]) +
       "\neach need to have a number and a sequence\n")

    try:
      pdb_file_start=int(alignment_lines[3].split()[0])
    except Exception:
      raise Sorry("\nSorry, the pdb_seq line %s " %(alignment_lines[3]) +
       " from your alignment file needs to have a number and a sequence")

    try:
      pdb_file_start=int(alignment_lines[4].split()[0])
    except Exception:
      raise Sorry("\nSorry, the pdb_seq line %s " %(alignment_lines[4]) +
       " from your alignment file needs to have a number and a sequence")

    pdb_file_alignment=alignment_lines[4].split()[1]
    seq_to_match_from_pdb,new_alignment_line=self.get_best_pdb_residue_list(
       pdb_file_alignment,pdb_file_start,out=out)
    if new_alignment_line: # we are going to use this later for editing
          pdb_file_alignment=new_alignment_line
          alignment_lines[4]=alignment_lines[4].split()[0]+" "+\
             new_alignment_line
          print("New alignment line 4 for PDB: \n%s" %(alignment_lines[4]), file=out)
          print("\n", file=out)

    # check the file, then really read it
    ok=True
    if len(alignment_lines) == 5 and alignment_lines[4].find("--") != 0: # just add it
      alignment_lines.append("--")

    if len(alignment_lines) != 6:
      print("Alignment file should have 6 lines. "+\
          "The file %s has these %s:\n\n%s" %(
        alignment_file,str(len(alignment_lines)),str("".join(alignment_lines))), file=out)
      ok=False

    if len(alignment_lines[0].split()) != 3 or \
          alignment_lines[0].split()[0] != "##":
      print("Alignment file line 1 should have 3 items: '##' 'xxxx' 'yyyy'", file=out)
      print("The file %s has %s: \n\n%s" %(alignment_file,
        str(len(alignment_lines[0].split())),alignment_lines[0]), file=out)
      ok=False

    for i in [3,4]:
      if len(alignment_lines[i].split()) != 2:
        print("Alignment file line %d should " %(i+1) + \
          "be a number and a string of characters", file=out)
        print("The file %s has %s: \n\n%s" %(alignment_file,
        str(len(alignment_lines[i].split())),alignment_lines[i]), file=out)
        ok=False

    try:
      seq_file_sequence=alignment_lines[3].split()[1]
      pdb_file_sequence=alignment_lines[4].split()[1]
    except Exception:
      ok=False
      raise Sorry("\nSorry,PDB or sequence-file sequences are missing?")

    # 2015-01-14:  put in check...
    # cannot handle situation where sequence file has
    # deletions relative to pdb file
    if seq_to_match_from_pdb.find(pdb_file_sequence.replace("-",""))<0:
      if getattr(self,
        'model_already_aligned_set_and_alignment_files_present',None):
        raise Sorry("Sorry, the parameter model_already_aligned was set to "+
        "\nTrue so your alignment file was ignored. "+
        "\n\nPlease set model_already_aligned=False so that your alignment "+
        "file can be read.\n")
      else:  #usual
        raise Sorry("Sorry, please supply an alignment file for this job.\n"+
         "Also specify force_alignment.\n"+
         "mr_rosetta cannot automatically generate alignments where there\n"+
         "are deletions required from the PDB file")

    if len(seq_file_sequence) != len(pdb_file_sequence):
      print("Length of sequences in alignment for seq_file (%d) and " %(
          len(seq_file_sequence)) + " PDB file (%d) do not match " %(
         len(pdb_file_sequence)), file=out)
      ok=False
    if not ok:
      raise Sorry ("Sorry, alignment file (\n%s\n)" %(
        "".join(alignment_lines)) +
        "should look something like this one:\n%s" %(standard_align) )

    # We expect line 4 of alignment file will be identical to seq from
    #  seq file...except possibly the start number is wrong.

    spl,start_pos,matched_sequence=self.get_matched_sequence(alignment_lines)
    print("Start position of target sequence in sequence file: ",start_pos,"\n", file=out)
    spl[0]=str(start_pos)
    line4=" ".join(spl)
    alignment_lines[3]=line4

    return alignment_lines,matched_sequence

  def get_matched_sequence(self,alignment_lines):
      line4=alignment_lines[3]
      spl=line4.split()
      if len(spl)!=2:
        raise Sorry("\nSorry, the sequence line %s " %(line4) +
       " from your alignment file needs to have a number and a sequence")
      try:
        n=int(spl[0])
      except Exception:
        raise Sorry("\nSorry, the sequence line %s " %(line4) +
       " from your alignment file needs to have a number and a sequence")
      alignment_seq=spl[1].replace("-","")
      # this must be PART of some member of self.sequence after
      #    removal of any "-"
      start_pos=-1
      matched_sequence=None
      for sequence in self.sequence:
        start_pos=sequence.find(alignment_seq)
        if start_pos>-1:
          matched_sequence=sequence
          break
      return spl,start_pos,matched_sequence

  def adjust_alignment_to_match_pdb(self,pdb_file_alignment,
        seq_to_match_from_pdb,out=sys.stdout,printout=True):
    # adjust pdb_file_alignment by editing to match seq_to_match_from_pdb, but
    # without changing the number of characters (i.e., only removing residues
    # that are not in seq_to_match_from_pdb and replacing with "-"

    # align seq_2 to seq_1, then replace seq_1 chars of X with ones from seq_2
    matches=0
    seq_1=pdb_file_alignment # fixed
    seq_2=seq_to_match_from_pdb # to be adjusted
    new_aa_for_fixed,score,break_string=\
      self.get_alignment(seq_1,seq_2,alignment_style='local',
          similarity_function="identity")
    # len(new_aa_for_fixed) is same as len(seq_1)
    if printout:
      print("\noriginal pdb_file_alignment:",pdb_file_alignment, file=out)
      print("new alignment              :",new_aa_for_fixed, file=out)
      if pdb_file_alignment.replace("-","") != new_aa_for_fixed.replace("-",""):
        print("NOTE: alignment has been edited", file=out)
      print(file=out)
    return new_aa_for_fixed,score

  def get_map_cc(self,params,s1=None,s2=None,map1=None,map2=None,labin_1=None,
      labin_2=None,out=sys.stdout):
    # get map correlation with phenix.get_cc_mtz_mtz between solutions s1 and s2
    if (s1 is not None and (
       not hasattr(s1,'map_coeffs') or not s1.map_coeffs)):
      return None
    if (s2 is not None and (
       not hasattr(s2,'map_coeffs') or not s2.map_coeffs)):
      return None
    if s1 and s2 and s1.get_space_group_symbol() != s2.get_space_group_symbol():
      return None
    if map1 is None and s1 is not None: map1=s1.map_coeffs
    if map2 is None and s2 is not None: map2=s2.map_coeffs
    if labin_1 is None and s1 is not None: labin_1=s1.labin_map_coeffs
    if labin_2 is None and s2 is not None: labin_2=s2.labin_map_coeffs
    if params.control.verbose:
      print("MAPS: ",map1,map2, file=out)
    from phenix.command_line.get_cc_mtz_mtz import get_cc_mtz_mtz
    quiet=not params.control.verbose
    from six.moves import cStringIO as StringIO
    f=StringIO()
    args=[map1,map2,
        "temp_dir="+str(params.directories.temp_dir),
        "output_dir="+str(params.directories.temp_dir)]
    if labin_1: args.append("labin_1="+labin_1)
    if labin_2: args.append("labin_2="+labin_2)
    g=get_cc_mtz_mtz(args,quiet=quiet,out=f)
    if params.control.verbose:
      print("CC VALUE: %6.2f " %(g.final_cc), file=out)
    return g.final_cc

  def read_input_solutions(self,params,out=sys.stdout,display=False):
      print("\n"+80*"=", file=out)
      print(" LOADING EXISTING SOLUTIONS", file=out)
      print(80*"=", file=out)
      print("Loading solutions from %s " %(params.input_files.mr_rosetta_solutions), file=out)
      # read results as pkl
      from libtbx import easy_pickle
      from phenix.rosetta.mr_rosetta import mr_rosetta_solution
      input_results=easy_pickle.load(params.input_files.mr_rosetta_solutions)
      print("RESULTS: ",input_results)
      if input_results is not None:
        if type(input_results) != type ([1,2,3]): input_results=[input_results]
        if hasattr(params.input_files,'ids_to_load') and \
           params.input_files.ids_to_load:
          print("Loading only the following solutions: %s" %(
            " ".join(list(params.input_files.ids_to_load))), file=out)
          self.input_solutions=[]
          for solution in input_results:
            if str(solution.id) in list(params.input_files.ids_to_load):
              self.input_solutions.append(solution)
        else:
          self.input_solutions=input_results
        for s in self.input_solutions:
            s.make_files_full_path()
        file=os.path.join(params.directories.temp_dir,'solutions_loaded.log')
        f=open(file,'w')
        print("Solutions loaded: ", file=f)
        for result in self.input_solutions:
          print(file=f)
          result.show_all(out=f)
          result.trace_parents(out=f,solutions=self.input_solutions)
          if hasattr(self,'highest_id'):
            self.highest_id=self.max(self.highest_id,result.id)
          else:
            self.highest_id=result.id
        f.close()
        print("\nLoaded %d previous solutions: \n(list is in %s) " %(
          len(self.input_solutions),file), file=out)
        if display:
           result_file,csv_file=self.set_result_file(params,out=out)
           print(open(file).read(), file=out)
           self.write_csv(params,saved_csv_file=csv_file,out=out)
        if params.non_user_params.write_local_files:
           self.results=self.input_solutions
           self.write_results(params,out=out,stage='stripped_results')

        # Extract self.model_info_object information
        if len(self.input_solutions)>=1:
          for s in self.input_solutions:
            if hasattr(s,'model_info_object') and \
                s.model_info_object is not None:
              self.extract_model_info(s,out=out)
              break
      else:
        print("No solutions found in %s" \
           %(params.input_files.mr_rosetta_solutions), file=out)

  def extract_model_info(self,solution,out=sys.stdout):
    if not solution.model_info_object: return
    self.model_info_object=solution.model_info_object
    for x in self.model_info_object.headers:
      if getattr(self.model_info_object,x) is not None:
        setattr(self,x,getattr(self.model_info_object,x))
    if solution.crystal_symmetry:  # 2011-01-07 use solution.crystal_symmetry if avail
      self.crystal_symmetry=solution.crystal_symmetry
    if hasattr(self,'crystal_symmetry') and self.crystal_symmetry:
      print("SET CRYSTAL SYMMETRY FROM INPUT SOLUTION: ",\
        iotbx.pdb.format_cryst1_record(crystal_symmetry=self.crystal_symmetry), file=out)
    if not hasattr(self,'ncs_object') or not self.ncs_object:
      self.ncs_object=solution.ncs_object

  def write_csv(self,params,saved_csv_file=None,out=sys.stdout):
    f=open(saved_csv_file,'w')
    print("\nWriting solutions as csv to %s " %(saved_csv_file), file=out)
    from phenix.rosetta.mr_rosetta import mr_rosetta_solution
    rs=mr_rosetta_solution(stage="mr_solution") # just a dummy
    possible_headers=rs.headers
    first_headers=["id","stage","r","r_free","best_mr_rescore_llg","rosetta_score"]
    headers=[]
    for h in first_headers+possible_headers:
      if not h in possible_headers: continue
      if h in headers: continue
      headers.append(h)
    print(",".join(headers), file=f)
    if not hasattr(self,'input_solutions'): self.input_solutions=[]
    for solution in self.input_solutions:
      value_list=[]
      for h in headers:
        if hasattr(solution,h):
           value=str(getattr(solution,h)).replace("\n"," ").replace(","," ")
        else:
           value='NA'
        value_list.append(value)
      print(",".join(value_list), file=f)
    f.close()

  def print_best_model_info(self,id,best_model,out=sys.stdout):
    print("ID: %d" %(id),": %s  %s" %(best_model.model,
       best_model.map_coeffs_file), file=out)
    if best_model.overall_model_cc:
      print(" Map-model CC: %8.2f" %(
        best_model.overall_model_cc), file=out)
    if best_model.r and best_model.r_free:
      print(" R/Rfree: %6.2f/%6.2f "% (
        best_model.r,best_model.r_free), file=out)

  def add_is_sub_process(self,edited_params): # mark as sub-process
    edited_params.non_user_params.is_sub_process=True
    # also turn off group_run_command if this  is a condor job
    if edited_params.control.condor or edited_params.control.one_subprocess_level:
       edited_params.control.group_run_command= \
         edited_params.control.single_run_command
       edited_params.control.condor=False
       edited_params.control.queue_commands=[] # 2011-07-01 was None...fixed
       edited_params.control.background=False

  def remove_stopwizard(self,params):
    file=os.path.join(self.params.directories.top_output_dir,'STOPWIZARD')
    if os.path.isfile(file):
      from phenix.autosol.delete_file import delete_file
      delete_file(file)

  def check_run_commands(self,params,out=sys.stdout):
    from phenix.autosol.check_run_command import check_run_command
    print("\nChecking single_run_command ...%s " %(self.single_run_command), end=' ', file=out)
    check_run_command(run_command=self.single_run_command,
       temp_dir=params.directories.temp_dir).run()
    print("OK", file=out)
    if self.group_run_command != self.single_run_command:
      print("\nChecking group_run_command ...%s " %(self.group_run_command), end=' ', file=out)
      check_run_command(run_command=self.group_run_command,
         temp_dir=params.directories.temp_dir).run()
      print("OK", file=out)

  def get_ncs_from_mr_models(self,params,mr_model_list,out=sys.stdout):
    # run simple_ncs_from_pdb to get ncs object..
    ncs_object_list=[]
    for model in mr_model_list:
      if model is None or not os.path.isfile(model):
        ncs_object_list.append(None)
      else:
        try:
          ncs_copies,ncs_object=self.get_ncs_from_model(params,
           model=model,out=out)
          if ncs_copies>1:
            ncs_object_list.append(ncs_object)
          else:
            ncs_object_list.append(None)
        except KeyboardInterrupt:
          self.stop_after_KeyboardInterrupt()
        except Exception:
          ncs_object_list.append(None)

    return ncs_object_list

  def get_ncs_from_model(self,params,model=None,out=sys.stdout):
    from phenix.command_line import simple_ncs_from_pdb
    file=os.path.join(params.directories.temp_dir,os.path.split(model)[-1][:-4]+"_ncs.log")
    f=open(file,'w')
    args_use=['pdb_in='+str(model)]
    args_use.append('temp_dir='+str(params.directories.temp_dir))
    args_use.append('verbose='+str(params.control.verbose))
    ncs_from_pdb=simple_ncs_from_pdb.run(args_use, log=f)
    ncs_object=ncs_from_pdb.get_ncs_object()
    if ncs_object.ncs_groups():
      ncs_copies=ncs_object.ncs_groups()[0].n_ncs_oper()
      print("\nNCS with %d copies found for %s\nLog file is %s\n " %(
        ncs_copies,model,file), file=out)
    else:
      ncs_copies=1
      print("\nNo NCS found for %s \nLog file is %s\n  " %(model,file), file=out)
    f.close()
    return ncs_copies,ncs_object

  def edit_model_for_repeat(self,params,temp=None,solution=None,out=sys.stdout):

    #Choose this model in 2 ways:
    #  1. take solution.model (an autobuild model) and use placed residues,
    #    refine this partial model
    #  2. take parent rosetta model, rs_refine, refine, and use that.
    #  Take the one with the better R-factor

    model_for_repeat=None
    type_of_model_for_repeat=None
    best_r=None
    autobuild_model=self.edit_autobuild_model(params,temp=temp,
      model=solution.model,out=out)
    if params.non_user_params.dummy_refinement: return autobuild_model

    if autobuild_model is not None:
      refined_prefix=autobuild_model[:-4]+"_ref"
      refined_autobuild_model,dummy_map_coeffs=self.run_refine(
         autobuild_model,refined_prefix,out=out,
         crystal_symmetry=solution.crystal_symmetry,
         refinement_params=params.input_files.refinement_params,
         correct_special_position_tolerance=
           params.non_user_params.correct_special_position_tolerance,
         skip_clash_guard=params.non_user_params.skip_clash_guard ,
         remove_clashing_residues=\
           params.refine_top_models.remove_clashing_residues,
         clash_cutoff=params.refine_top_models.clash_cutoff,
         Facts={}) # 2013-09-12 # 2016-02-27
      if refined_autobuild_model is not None and (
           best_r is None or self.r<best_r):
        model_for_repeat=refined_autobuild_model
        type_of_model_for_repeat="refined_autobuild_model"
        best_r=self.r
    else:
      print("AutoBuild model not suitable for repeat (no placed residues)\n", file=out)


    rosetta_model=solution.parent_model
    if params.non_user_params.dummy_refinement: return rosetta_model
    if rosetta_model is not None:
      print("\nRefining rosetta model using map from autobuild", file=out)
      # First rs_refine against autobuild map (solution.map_coeffs)

      # Remove free set:
      new_mtz=solution.map_coeffs[:-4]+"_nf.mtz"
      new_labin=self.remove_free(mtz=solution.map_coeffs,
          labin=solution.labin_map_coeffs,
          free=params.input_files.data,
          labin_free=params.input_files.labin,
          new_mtz=new_mtz,
          test_flag_value=params.control.test_flag_value,
          out=out,temp_dir=params.directories.temp_dir)
      self.wait_for_file_to_appear(new_mtz,
       max_ticks=params.control.max_wait_time,allow_failure=False)

      # Real-space refine:
      new_pdb=rosetta_model[:-4]+"_rs.pdb"
      print("\nRunning real-space refinement on %s to yield %s  (score is CC)"%(
        rosetta_model,new_pdb), file=out)
      from phenix.utilities.rs_refine import rs_refine
      rs=rs_refine(setup=True,mtz_in=new_mtz,labin=new_labin)
      rs.refine_pdb(pdb=rosetta_model,pdb_out=new_pdb,
        use_space_group_from_pdb=True, # 2013-08-07
       )


      # Now refine in standard way
      refined_prefix=new_pdb[:-4]+"_ref"
      refined_rosetta_model,dummy_map_coeffs=self.run_refine(
         new_pdb,refined_prefix,out=out,
         crystal_symmetry=solution.crystal_symmetry,
         refinement_params=params.input_files.refinement_params,
         correct_special_position_tolerance=
           params.non_user_params.correct_special_position_tolerance,
         skip_clash_guard=params.non_user_params.skip_clash_guard ,
         remove_clashing_residues=\
           params.refine_top_models.remove_clashing_residues,
         clash_cutoff=params.refine_top_models.clash_cutoff,
         Facts={}) # 2013-09-12 # 2016-02-27
      if refined_rosetta_model is not None and (
           best_r is None or self.r<best_r):
        model_for_repeat=refined_rosetta_model
        type_of_model_for_repeat="refined_rosetta_model"
        best_r=self.r
    if model_for_repeat:
      print("\n Final model for repeat is %s: %s with R=%6.2f" %(
        type_of_model_for_repeat,model_for_repeat,best_r), file=out)
    return model_for_repeat

  def edit_autobuild_model(self,params,temp=None,model=None,
      out=sys.stdout):
    from phenix.autosol.atomname import atomname
    from phenix.autosol.trim_file_name import trim_file_name
    model_name=trim_file_name(model).trimmed_file

    output_pdb=os.path.join(temp,model_name[:-4]+"_tr.pdb")
    pdb_input = get_pdb_inp(file_name=model)
    hierarchy = pdb_input.construct_hierarchy()
    # take all residues up to number in seq file (placed residues only)
    max_seq=1
    for x in self.sequence_by_group:
      if len(x)>max_seq: max_seq=len(x)
    print("\nEditing %s keeping only placed residues (resno < %d) " %(
      model,max_seq), file=out)
    f=open(output_pdb,'w')
    self.get_crystal_symmetry_from_pdb(model=model,pdb_input=pdb_input,out=out)
    print(iotbx.pdb.format_cryst1_record(
         crystal_symmetry=self.crystal_symmetry), file=f)
    residues_kept=0
    residues_skipped=0
    for model in hierarchy.models()[:1]:
        for chain in model.chains():
          for residue_group in chain.residue_groups():
            if residue_group.resseq_as_int() > max_seq:
               residues_skipped+=1
               continue
            is_water=False
            for r in residue_group.unique_resnames():
              if atomname().is_water(r):
                is_water=True
            if is_water:
               continue
            if residue_group.resseq_as_int() > max_seq:
               residues_skipped+=1
               continue
            residues_kept+=1
            for atom_group in residue_group.atom_groups():
              for atom in atom_group.atoms():
                print(atom.format_atom_record(), file=f)
    f.close()
    print("Total of %d residues kept and %d skipped; written to %s "%(
       residues_kept,residues_skipped,output_pdb))
    if residues_kept<1:
      return None
    else:
      return output_pdb

  def replace_two_sep(self,path,trailing_sep=True):
    sep=os.path.sep
    two_sep=sep+sep
    edited_path=path
    while two_sep in edited_path:
      edited_path=edited_path.replace(two_sep,sep)
    if trailing_sep and edited_path[-1]!=sep:  # add trailing separator
      edited_path=edited_path+sep
    return edited_path

  def make_relative(self,path,local_path):
    # remove local path from path if possible
    edited_path=self.replace_two_sep(path,trailing_sep=False)
    edited_local_path=self.replace_two_sep(local_path)
    if edited_path.find(edited_local_path)==0:
      edited_path=edited_path.replace(edited_local_path,"")
      return edited_path
    return path

  def strings_from_anything(self,anything_list):
    s_list=[]
    for f in anything_list:
      s=str(f)
      s_list.append(s)
    return " ".join(s_list)

  def get_copies_to_find(self,params,
    ncs_copies=0,model=None,n_fixed=0,out=sys.stdout):
    print("\nGuessing number of copies to find\n", file=out)
    copies_in_search_model,dummy_ncs=\
        self.get_ncs_from_model(params,model=model,out=out)
    copies_of_search_model=n_fixed
    copies_total=ncs_copies
    copies_to_find=copies_total-copies_of_search_model*copies_in_search_model
    copies_of_search_model_to_find=copies_to_find//copies_in_search_model
    print("Copies in search model: %d" %(copies_in_search_model), file=out)
    print("Copies of search model already placed: %d" %(copies_of_search_model), file=out)
    print("Total copies in asymmetric unit: %d" %(copies_total), file=out)
    print("Copies to find: %d" %(copies_to_find), file=out)
    print("Copies of search model to find: %d\n" %(copies_of_search_model_to_find), file=out)
    if copies_of_search_model_to_find*copies_in_search_model != copies_to_find:
      print("\nWARNING: copies needed is not a multiple of copies in"+\
        " search model\n", file=out)
    if copies_of_search_model_to_find==0:
      raise Sorry("Number of copies to look for is 0...something must not be set up quite correctly.")
    return copies_of_search_model_to_find

  def check_nproc_and_background(self,params,out=sys.stdout):
    if params.control.background is None and params.control.nproc > 1:
      if params.control.group_run_command.replace(" ","").lower()=="sh":
        print("Running in background as nproc=%d and group_run_command=%s" %(
          params.control.nproc,params.control.group_run_command))
        params.control.background=True

  def set_up_weights(self,params,out=sys.stdout):
    # set up weights file for rosetta if needed
    if params.rosetta_modeling.weights_file:  # use weights provided
      if params.rosetta_modeling.include_solvation_energy==False:
        raise Sorry("Sorry, you cannot specify a weights_file and "+
           "include_solvation_energy=False.\nPlease either supply an edited "+
           "weights_file or just turn off include_solvent_energy.\n" )
      if not os.path.isfile(params.rosetta_modeling.weights_file):
        raise Sorry("Sorry, the weights_file %s is missing?" %(
          params.rosetta_modeling.weights_file))
      else: # all ok
        print("\nUsing the file %s for Rosetta weights\n" %(
          params.rosetta_modeling.weights_file), file=out)
    elif params.rosetta_modeling.include_solvation_energy==False: # set weights
      orig_weights_file=os.path.join(self.rosetta_database_dir,'scoring',
         'weights','score12_full.wts')
      if not os.path.isfile(orig_weights_file):
        raise Sorry("Sorry the weights file %s is missing?") %(
          orig_weights_file)
      new_weights_file=os.path.join(params.directories.temp_dir,
         'rosetta_weights.wts')
      f=open(new_weights_file,'w')
      for line in open(orig_weights_file).readlines():
        if line.find('fa_sol ')>-1:
          line="fa_sol 0.0"
        print(line.rstrip(), file=f)
      f.close()
      self.wait_for_file_to_appear(new_weights_file,
       max_ticks=params.control.max_wait_time,allow_failure=False)
      params.rosetta_modeling.weights_file=new_weights_file
      params.rosetta_modeling.include_solvation_energy=None
      print("\nEditing the weights file %s" %(orig_weights_file) +\
        "\nto remove solvation energy. New weights: %s" %(
          params.rosetta_modeling.weights_file), file=out)

  def set_space_group_object(self,sg_symbol):
    # return None if not present, otherwise sg object
    if sg_symbol:
      from cctbx import sgtbx
      return sgtbx.space_group_info(sg_symbol)
    else:
      return  None

  def get_ncs_copies_from_params(self,params):
    # allows it to be params.crystal_info.ncs_copies or
    # params.crystal_info.ncs_copies[0]
    if not hasattr(params,'crystal_info'): return
    if not hasattr(params.crystal_info,'ncs_copies'): return
    if not params.crystal_info.ncs_copies: return

    if type (params.crystal_info.ncs_copies)==type([1,2,3]):
      ncs_copies=params.crystal_info.ncs_copies[0]
    else:
      ncs_copies=params.crystal_info.ncs_copies
    return ncs_copies

  def set_initial_autobuild_args(self,params,solution=None,out=sys.stdout,
      data=None):
    if data is None and solution is not None:
       data=solution.data
    if data is None:
      data=params.input_files.data
    if data is None:
      raise Sorry("Sorry, need data for set_initial_autobuild_args")
    args=["data="+str(data),
          "top_output_dir="+str(params.directories.top_output_dir),
          "nproc="+str(params.control.nproc),
          "run_command="+str(params.control.group_run_command),
          "background="+str(params.control.background),
          "ignore_errors_in_subprocess="+
             str(params.control.ignore_errors_in_subprocess)
         ]
    if hasattr(params.input_files,'seq_file'):
       args.append("seq_file="+str(params.input_files.seq_file))

    if hasattr(params.crystal_info,'space_group') and \
        params.crystal_info.space_group:
       args.append("space_group="+str(params.crystal_info.space_group))
    if params.crystal_info.resolution:
      args.append("resolution="+str(params.crystal_info.resolution))

    if hasattr(params.crystal_info,'solvent_fraction') and \
         params.crystal_info.solvent_fraction:
      params.crystal_info.solvent_fraction=\
       self.put_solvent_in_range(params.crystal_info.solvent_fraction,out=out)

      args.append("solvent_fraction="+str(params.crystal_info.solvent_fraction))

    if self.params.input_files.labin:
      labels=self.get_labels_from_labin(
          self.params.input_files.labin,out=out)
      if labels is not None:
        args.append("input_labels="+labels)

    if params.control.verbose:
      args.append("verbose=true")
      print("Keywords: %s" %(str(args)), file=out)
    if params.control.debug:
      args.append("debug=true")

    ncs_copies=self.get_ncs_copies_from_params(params)
    print("NCS copies obtained from parameters:",ncs_copies)

    if ncs_copies:
      args.append("ncs_copies="+str(ncs_copies))
    elif solution and solution.ncs_copies: # use if available
      args.append("ncs_copies="+str(solution.ncs_copies))

    if hasattr(params.non_user_params,'correct_special_position_tolerance'):
      if params.non_user_params.correct_special_position_tolerance is not None:
       args.append("correct_special_position_tolerance="+float(
         params.non_user_params.correct_special_position_tolerance))
      if params.non_user_params.skip_clash_guard is not None:
       args.append("skip_clash_guard=True")
    if hasattr(params,'refinement') and \
        hasattr(params.refinement,'correct_special_position_tolerance'):
      if params.refinement.correct_special_position_tolerance is not None:
       args.append("correct_special_position_tolerance="+float(
         params.refinement.correct_special_position_tolerance))
      if params.refinement.skip_clash_guard is not None:
       args.append("skip_clash_guard=True")

    if hasattr(params.crystal_info,'chain_type'):
       args.append("chain_type=%s" %(params.crystal_info.chain_type))

    if hasattr(params.control,'add_double_quotes_in_condor'):
       args.append("add_double_quotes_in_condor=%s" %(
          params.control.add_double_quotes_in_condor))

    if hasattr(params.control,'condor'):
       args.append("condor=%s" %(params.control.condor))

    if hasattr(params.control,'condor_universe'):
       args.append("condor_universe=%s" %(params.control.condor_universe))

    if hasattr(params.control,'queue_commands'):
      for x in params.control.queue_commands:
        args.append("queue_commands=%s" %(x))


    return args
