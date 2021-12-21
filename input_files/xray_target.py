
from __future__ import division
from mmtbx.maps import composite_omit_map
import mmtbx.command_line
import mmtbx.maps.fem
import mmtbx.utils
import iotbx.map_tools
from cctbx import maptbx
from cctbx import xray
from scitbx.array_family import flex
from libtbx.str_utils import make_header
from libtbx.utils import Usage, null_out
from libtbx import Auto
import time
import os
import sys

map_type_phil_str = """
map_type = *2mFo-DFc mFo dm prime_and_switch feature_enhanced composite_omit
  .type = choice
  .short_caption = Real-space refinement map type
  .help = This controls the procedure for generating a map for Rosetta \
    real-space refinement.  All map types other than 'mFo' will essentially \
    be similar to 2mFo-DFc, but with different processing steps.
"""

def master_phil () :
  return mmtbx.command_line.generate_master_phil_with_inputs(
    enable_automatic_twin_detection=False,
    phil_string="""
output {
  write_grads = True
    .type = bool
  write_eff = True
    .type = bool
  prefix = gradients
    .type = str
  serial = 1
    .type = int
  output_dir = None
    .type = path
  ccp4_map {
    write_map = False
      .type = bool
    %s
    map_file_name = phenix_map.ccp4
      .type = path
    exclude_free_r_flags = True
      .type = bool
    fill_missing_f_obs = False
      .type = bool
    grid_resolution_factor = 1/3
      .type = float
    sharpen = True
      .type = bool
    sharpen_b = 0
      .type = float
    composite_omit = False
      .type = bool
  }
}
options {
  bulk_solvent_and_scale = True
    .type = bool
  optimize_mask = False
    .type = bool
  force_bss = False
    .type = bool
  remove_outliers = True
    .type = bool
  target_name = *ml lsq mlhl
    .type = choice(multi=False)
  twin_law = None
    .type = str
  compute_gradients = True
    .type = bool
}
resolve {
  density_modify = False
    .type = bool
  prime_and_switch = False
    .type = bool
  solvent_content = 0.5
    .type = float
  mask_cycles = 5
    .type = int
  minor_cycles = 10
    .type = float
}
include scope mmtbx.command_line.fmodel.fmodel_from_xray_structure_params
bss {
  include scope mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params
}
""" % map_type_phil_str)

class xray_target (object) :
  def __init__ (self,
                params,
                xray_structure,
                pdb_hierarchy, # used for NCS averaging later
                f_obs,
                r_free_flags,
                hl_coeffs=None,
                out=sys.stdout,
                ignore_bulk_solvent=False) :
    self.params = params
    self.out = out
    t1 = time.time()
    target_name = params.options.target_name
    if (target_name == "lsq") :
      target_name = "ls_wunit_k1"
    if (hl_coeffs is not None) :
      hl_coeffs = hl_coeffs.common_set(f_obs)
      if (hl_coeffs.indices().size() != f_obs.indices().size()) :
        raise RuntimeError("Different reflection counts for data and "+
          "Hendrickson-Lattman coefficients!")
    fmodel = mmtbx.utils.fmodel_manager(
      f_obs=f_obs,
      r_free_flags=r_free_flags,
      xray_structure=xray_structure,
      sf_and_grads_accuracy_params=params.structure_factors_accuracy,
      mask_params=params.mask,
      target_name=target_name,
      hl_coeff=hl_coeffs,
      twin_law=params.options.twin_law)
    #fmodel_info = fmodel.info()
    #fmodel_info.show_rfactors_targets_scales_overall(out=out)
    if params.options.remove_outliers :
      print >> out, ""
      fmodel.remove_outliers(log=out)
    self.fmodel = fmodel
    if params.options.bulk_solvent_and_scale :
      self.fmodel.update_all_scales(
        params=self.params.bss,
        fast=True,
        optimize_mask=self.params.options.optimize_mask,
        show=True)
      fmodel_info = fmodel.info()
      fmodel_info.show_rfactors_targets_scales_overall(out=out)
    else :
      make_header("Update X-ray structure", out=out)
      #fmodel.apply_back_b_iso()
      fmodel.update_xray_structure(update_f_mask=True)
    self.xray_structure = fmodel.xray_structure
    self.flag_apply_shifts = False
    self._update_target_functor = False

  def optimize_mask (self) :
    self.fmodel.optimize_mask()

  #XXX we do not want to use fmodel.update_all_scales() here, because this
  # may change the size of f_obs.
  def optimize_mask_and_update_solvent_and_scale (self) :
    make_header("Bulk solvent and scaling", out=self.out)
    self.fmodel.update_all_scales(
      params=self.params.bss,
      remove_outliers=False, # see comment above
      fast=True)
    self.fmodel.optimize_mask(params=self.params.bss)
    #self._update_target_functor = True

  def update_solvent_and_scale (self):
    make_header("Bulk solvent and scaling", out=self.out)
    self.fmodel.update_all_scales(
      params=self.params.bss,
      remove_outliers=False, # see comment above
      fast=True)
    #self._update_target_functor = True

  def update_fmask (self) :
    self.fmodel.update_xray_structure(
      xray_structure=self.xray_structure,
      update_f_calc=False,
      update_f_mask=True)

  def prepare_for_minimization (self) :
    make_header("Preparing for minimization", out=self.out)
    xrs = self.xray_structure
    xrs.scatterers().flags_set_grads(state=False)
    selection = flex.bool(xrs.sites_cart().size(), True)
    #print selection.iselection().size()
    xrs.scatterers().flags_set_grad_site(iselection=selection.iselection())
    self.x = flex.double(xrs.n_parameters(), 0)
    self._scatterers_start = xrs.scatterers()
    self.fmodel.update_xray_structure(xray_structure=self.xray_structure)
    self.target_functor = self.fmodel.target_functor()
    self.target_functor.prepare_for_minimization()
    self.fmodel.info().show_targets(out=self.out, text="Target values")

  def compute_target (self, compute_gradients=True, use_vec3_array=False) :
    if (getattr(self, "target_functor", None) is None) :
      raise RuntimeError("Please call prepare_for_minimization() before "+
        "computing target and gradients.")
    if (self.flag_apply_shifts) :
      self.apply_shifts()
    out = self.out
    fmodel = self.fmodel
    fmodel.update_xray_structure(xray_structure=self.xray_structure,
      update_f_calc=True)
    #if (self._update_target_functor) :
    #  self.target_functor = self.fmodel.target_functor()
    #  self.target_functor.prepare_for_minimization()
    t1 = time.time()
    tfx_r = self.target_functor(compute_gradients=compute_gradients)
    self.f = tfx_r.target_work()
    if compute_gradients :
      grads = tfx_r.gradients_wrt_atomic_parameters().packed()
      assert (grads.size() == self.x.size())
      if (use_vec3_array) :
        self.g = flex.vec3_double(grads)
      else :
        self.g = grads
    else :
      if (use_vec3_array) :
        self.g = flex.vec3_double(int(self.x.size() / 3), (0.0,0.0,0.0))
      else :
        self.g = flex.double(self.x.size(), 0.0)
    t2 = time.time()
    return self.f, self.g

  def compute_functional_and_gradients (self) :
    self.apply_shifts()
    return self.compute_target()

  def compute_functional_and_gradients_rosetta (self) :
    f, g = self.compute_target(compute_gradients=True)
    return f, list(g)

  def target (self) :
    return self.f

  def gradients (self) :
    return self.g

  def r_work (self) :
    return self.fmodel.r_work()

  def r_free (self) :
    return self.fmodel.r_free()

  def apply_shifts (self) :
    xrs = self.xray_structure
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = xrs.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    site_symmetry_table = xrs.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      scatterers_shifted[i_seq].site = crystal.correct_special_position(
        crystal_symmetry = self.xray_structure,
        special_op       = site_symmetry_table.get(i_seq).special_op(),
        site_frac        = scatterers_shifted[i_seq].site,
        site_label       = scatterers_shifted[i_seq].label,
        tolerance        = self.correct_special_position_tolerance)
    xrs.replace_scatterers(scatterers = scatterers_shifted)
    return None

  def update_sites (self, sites_cart) :
    self.xray_structure.set_sites_cart(sites_cart)
    self.fmodel.update_xray_structure(
      xray_structure=self.xray_structure,
      update_f_calc=False,
      update_f_mask=True)

  def update_sites_1d (self, sites_cart_as_1d) :
    if isinstance(sites_cart_as_1d, list) :
      assert (len(sites_cart_as_1d) == self.x.size())
      sites_cart_as_1d = flex.double(sites_cart_as_1d)
    else :
      assert (sites_cart_as_1d.size() == self.x.size())
    sites_cart = flex.vec3_double(sites_cart_as_1d)
    self.update_sites(sites_cart)

  def set_shifts (self, site_shifts) :
    assert (site_shifts.size() == self.x.size())
    self.x = site_shifts

  def clean_up_after_minimization (self) :
    self.apply_shifts()
    info = self.fmodel.info()
    info.show_targets(out=self.out, text="Target values after minimization")

  def update_pdb_hierarchy (self, pdb_hierarchy) :
    pdb_hierarchy.atoms().set_xyz(self.xray_structure.sites_cart())

  def compute_ccp4_map (self,
                        map_type="2mFo-DFc",
                        file_name="2mFo-DFc.ccp4") :
    params = self.params
    fill_missing_f_obs = params.output.ccp4_map.fill_missing_f_obs
    grid_resolution_factor = params.output.ccp4_map.grid_resolution_factor
    sharpen = params.output.ccp4_map.sharpen
    sharpen_b = params.output.ccp4_map.sharpen_b
    self.fmodel.update_all_scales(
      params=self.params.bss,
      fast=True,
      optimize_mask=False, #self.params.options.optimize_mask,
      show=False)
    # XXX does the bulk solvent + scaling need to be updated again here?
    map_type_str = "2mFo-DFc"
    if (map_type == "mFo") :
      map_type_str = "mFo"
    if (map_type == "composite_omit") :
      crystal_gridding = self.fmodel.f_obs().crystal_gridding(
        resolution_factor=grid_resolution_factor,
        symmetry_flags=maptbx.use_space_group_symmetry,
        assert_shannon_sampling=True)
      map_coeffs = composite_omit_map.run(
        crystal_gridding=crystal_gridding,
        fmodel=self.fmodel,
        map_type=map_type,
        resolution_factor=grid_resolution_factor,
        exclude_free_r_reflections=True,
        fill_missing=fill_missing_f_obs,
        log=null_out()).map_coefficients
    elif (map_type == "feature_enhanced") :
      map_coeffs = mmtbx.maps.fem.run(
        fmodel=self.fmodel,
        signal_threshold=0.5,
        sharp=True).mc_result
      map_coeffs, r_free_flags = map_coeffs.common_sets(
        other=self.fmodel.r_free_flags())
      map_coeffs = map_coeffs.select(~(r_free_flags.data()))
    else :
      map_coeffs = self.fmodel.map_coefficients(
        map_type=map_type_str,
        exclude_free_r_reflections=True,
        fill_missing=fill_missing_f_obs)
    if (sharpen == True) or ((sharpen is Auto) and
                             (map_type != "feature_enhanced"))  :
      # XXX following the convention used in AutoBuild and related tools,
      # positive will sharpen, negative will blur - note that this is the
      # opposite of the CCTBX API (and possibly phenix.refine/phenix.maps)
      b_sharp = 10 * map_coeffs.d_min()
      if (sharpen_b != 0) and (sharpen_b is not None) :
        b_sharp = sharpen_b
      print >> self.out, "Sharpening with B=-%g" % b_sharp
      map_coeffs = map_coeffs.apply_debye_waller_factors(b_iso=-b_sharp)
    if map_coeffs.anomalous_flag() :
      map_coeffs = map_coeffs.average_bijvoet_mates()
    if ((map_type in ["dm", "prime_and_switch"]) and
        (self.fmodel.twin_law is None)) :
      print >> self.out, "Running RESOLVE..."
      from phenix.command_line import prime_and_switch_map
      prime_and_switch_map.run_resolve_dm(
        map_coeffs=map_coeffs,
        fmodel=self.fmodel,
        params=self.params.resolve,
        prime_and_switch=(map_type == "prime_and_switch"),
        exclude_free_r_reflections=exclude_free_r_flags,
        out=null_out(),
        verbose=False)
    fft_map = map_coeffs.fft_map(resolution_factor=grid_resolution_factor)
    fft_map.apply_sigma_scaling()
    iotbx.map_tools.write_ccp4_map(
      sites_cart=self.xray_structure.sites_cart(),
      unit_cell=fft_map.unit_cell(),
      map_data=fft_map.real_map(),
      n_real=fft_map.n_real(),
      file_name=file_name)

  def write_map (self, file_name=None) :
    params = self.params
    if (file_name is None) :
      file_name = params.output.ccp4_map.map_file_name
    self.compute_ccp4_map(
      map_type=params.output.ccp4_map.map_type,
      file_name=file_name)
    print >> self.out, "CCP4 map:"
    print >> self.out, "  %s" % params.output.ccp4_map.map_file_name

  def write_files (self) :
    out = self.out
    params = self.params
    fmodel = self.fmodel
    print >> out, ""
    if (params.output.output_dir is None) :
      params.output.output_dir = os.getcwd()
    file_base = os.path.join(params.output.output_dir,
      "%s_%d" % (params.output.prefix, params.output.serial))
    params.output.serial += 1
    if (not params.options.force_bss) :
      params.options.bulk_solvent_and_scale = False
    if params.output.write_eff :
      f = open("%s.eff" % file_base, "w")
      final_phil = master_phil().format(python_object=params)
      final_phil.show(out=f)
      f.close()
      print >> out, "Parameters for next run:"
      print >> out, "  %s.eff" % file_base
    if params.output.write_grads :
      grads_out = open("%s.xyz" % file_base, "w")
      print >> grads_out, "# TARGET = %.12g" % self.f
      if isinstance(self.g, flex.double) :
        grads = flex.vec3_double(self.g)
      else :
        grads = self.g
        for site_grads in grads :
          print >> grads_out, "%e %e %e" % site_grads
      grads_out.close()
      print >> out, "Gradients:"
      print >> out, "  %s.xyz" % file_base
    if params.output.ccp4_map.write_map :
      self.write_map()

def run (args, out=sys.stdout) :
  if (len(args) == 0) :
    raise Usage("""
phenix.compute_xray_gradients [model.pdb] [data.mtz] [params.eff] [options...]

This program will write out a parameter (.eff) file when it is finished,
which can be used for the next round without bulk-solvent correction.
""")
  elif (args == ["--options"]) or (args == ["--help"]) :
    master_phil().show()
    return None
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=False,
    create_fmodel=False)
  params = cmdline.params
  f_obs = cmdline.f_obs
  r_free_flags = cmdline.r_free_flags
  xray_structure = cmdline.xray_structure
  hl_coeffs = None
  if (params.options.target_name == "mlhl") :
    for array in cmdline.miller_arrays :
      if (array.is_hendrickson_lattman_array()) :
        hl_coeffs = array
        break
    if (hl_coeffs is None) :
      raise Sorry("Hendrickson-Lattman coefficients required when "+
        "target_name=mlhl.")
  target_evaluator = xray_target(
    params=params,
    xray_structure=cmdline.xray_structure,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    f_obs=cmdline.f_obs,
    r_free_flags=cmdline.r_free_flags,
    hl_coeffs=hl_coeffs,
    out=out)
  target_evaluator.prepare_for_minimization()
  target_evaluator.compute_target(
    compute_gradients=params.options.compute_gradients,
    use_vec3_array=True)
  target_evaluator.write_files()
  return target_evaluator

if __name__ == "__main__" :
  run(sys.argv[1:])
