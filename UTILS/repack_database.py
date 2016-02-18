#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Repacking Instaseis databases.

Requires click, h5py, and numpy.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2016
    Simon Stähler (staehler@geophysik.uni-muenchen.de), 2016

:license:
    GNU Lesser General Public License, Version 3 [non-commercial/academic use]
    (http://www.gnu.org/copyleft/lgpl.html)
"""
import os

import click
import numpy as np


def maybe_encode(string, encoding='ascii'):
    try:
        return string.encode(encoding)
    except AttributeError:
        return string
    except UnicodeEncodeError:
        return string

def unroll_and_merge_netcdf4(filenames, output_folder):
    """
    Completely unroll and merge both files.
    """

    import netCDF4
    # Find MZZ, MXX_P_MYY, MXZ_MYZ, MXY_MXX_M_MYY directories
    if len(filenames) == 4:
        filenames = [os.path.normpath(_i) for _i in filenames]
        mzz = [_i for _i in filenames if "MZZ" in _i]
        mxx = [_i for _i in filenames if "MXX_P_MYY" in _i]
        mxz = [_i for _i in filenames if "MXZ_MYZ" in _i]
        mxy = [_i for _i in filenames if "MXY_MXX_M_MYY" in _i]
        assert len(mzz) == 1
        assert len(mxx) == 1
        assert len(mxz) == 1
        assert len(mxy) == 1
        mzz = mzz[0]
        mxx = mxx[0]
        mxz = mxz[0]
        mxy = mxy[0]
        assert os.path.exists(mzz)
        assert os.path.exists(mxx)
        assert os.path.exists(mxz)
        assert os.path.exists(mxy)
        f_in_1 = netCDF4.Dataset(mzz, 'r')
        f_in_2 = netCDF4.Dataset(mxx, 'r')
        f_in_3 = netCDF4.Dataset(mxz, 'r')
        f_in_4 = netCDF4.Dataset(mxy, 'r')

    elif len(filenames) == 2:
        pz = [_i for _i in filenames if "PZ" in _i]
        px = [_i for _i in filenames if "PX" in _i]
        assert len(pz) == 1
        assert len(px) == 1
        pz = pz[0]
        px = px[0]
        assert os.path.exists(pz)
        assert os.path.exists(px)
        f_in_1 = netCDF4.Dataset(pz, 'r')
        f_in_2 = netCDF4.Dataset(px, 'r')

    output_filename = os.path.join(output_folder, "merged_instaseis_db.nc4")
    assert not os.path.exists(output_filename)

    assert not os.path.exists(output_filename)

    try:
        f_out = netCDF4.Dataset(output_filename, 'w', format='NETCDF4')

        # Copy attributes from the vertical file.
        for name in f_in_1.ncattrs():
            #value = f_in_1.getncattr(name)
            value = getattr(f_in_1, name)
            print name, value
            #f_out.setncattr(name, value)
            setattr(f_out, name, maybe_encode(value))

        f_out.setncattr('nsim', len(filenames))

        for name, dimension in f_in_1.dimensions.iteritems():
            if not dimension.isunlimited():
                f_out.createDimension(name, len(dimension))
            else:
                f_out.createDimension(name, None)

        # Create Mesh group and copy mesh variables
        f_out.createGroup('Mesh')
        for name, dimension in f_in_1['Mesh'].dimensions.iteritems():
            if not dimension.isunlimited():
                f_out['Mesh'].createDimension(name, len(dimension))
            else:
                f_out['Mesh'].createDimension(name, None)

        for name, variable in f_in_1['Mesh'].variables.iteritems():
            f_out['Mesh'].createVariable(name, variable.datatype,
                                         variable.dimensions)
            f_out['Mesh'].variables[name][:] = f_in_1['Mesh'].variables[name][:]

        # Copy source time function variables from Surface group
        #f_out.createGroup('Surface')
        #for name, dimension in f_in_1['Surface'].dimensions.iteritems():
        #    if not dimension.isunlimited():
        #        f_out['Surface'].createDimension(name, len(dimension))
        #    else:
        #        f_out['Surface'].createDimension(name, None)

        for name, variable in f_in_1['Surface'].variables.iteritems():
            if name in ['stf_dump', 'stf_d_dump']:
                f_out.createVariable(name, variable.datatype,
                                     variable.dimensions)
                f_out.variables[name][:] = f_in_1['Surface'].variables[name][:]

        # Create a new array but this time in 5D. The first dimension
        # is the element number, the second and third are the GLL
        # points in both directions, the fourth is the time axis, and the
        # last the displacement axis.
        npts = f_in_1.getncattr("number of strain dumps")
        number_of_elements = f_in_1.getncattr("nelem_kwf_global")
        npol = f_in_1.getncattr("npol")

        # Get datasets and the dtype.
        if len(filenames) == 2:
            meshes = [
                f_in_1["Snapshots"]["disp_s"], # PZ
                f_in_1["Snapshots"]["disp_z"],
                f_in_2["Snapshots"]["disp_s"], # PX
                f_in_2["Snapshots"]["disp_p"],
                f_in_2["Snapshots"]["disp_z"]]
        elif len(filenames) == 4:
            meshes = [
                f_in_1["Snapshots"]["disp_s"],  # MZZ
                f_in_1["Snapshots"]["disp_z"],
                f_in_2["Snapshots"]["disp_s"],  # MXX + MYY
                f_in_2["Snapshots"]["disp_z"],
                f_in_3["Snapshots"]["disp_s"],  # MXZ / MYZ
                f_in_3["Snapshots"]["disp_p"],
                f_in_3["Snapshots"]["disp_z"],
                f_in_4["Snapshots"]["disp_s"],  # MXY / MXX - MYY
                f_in_4["Snapshots"]["disp_p"],
                f_in_4["Snapshots"]["disp_z"]]

        dtype = meshes[0].dtype

        nvars = len(meshes)

        dim_elements = f_out.createDimension('elements', number_of_elements)
        dim_ipol = f_out.createDimension('ipol', npol + 1)
        dim_jpol = f_out.createDimension('jpol', npol + 1)
        dim_nvars = f_out.createDimension('variables', nvars)
        dim_snaps = f_out.dimensions['snapshots']

        print f_out.dimensions

        ds_o = f_out.createVariable(varname="merged_snapshots",
                                    dimensions=(dim_elements.name,
                                                dim_nvars.name,
                                                dim_jpol.name,
                                                dim_ipol.name,
                                                dim_snaps.name),
                                    datatype=dtype, contiguous=True)

                                    # Old order (Instaseis):
                                    # dimensions=(dim_elements.name,
                                    #             dim_snaps.name,
                                    #             dim_ipol.name,
                                    #             dim_jpol.name,
                                    #             dim_nvars.name),
        print ds_o

        utemp = np.zeros((nvars, npol + 1, npol + 1, npts),
                         dtype=dtype)
        print utemp.shape

        # Now it becomes more interesting and very slow.
        sem_mesh = f_in_1["Mesh"]["sem_mesh"]
        with click.progressbar(range(number_of_elements),
                               length=number_of_elements,
                               label="\t  ") as gll_idxs:
            for gll_idx in gll_idxs:
                gll_point_ids = sem_mesh[gll_idx]

                # Load displacement from all GLL points.
                for ivar, var in enumerate(meshes):
                    # The list of ids we have is unique but not sorted.
                    ids = gll_point_ids.flatten()
                    s_ids = np.sort(ids)
                    temp = var[:, s_ids]
                    for ipol in range(npol + 1):
                        for jpol in range(npol + 1):
                            idx = ipol * (npol + 1) + jpol
                            # ipol, jpol, nsnap, nvar
                            # ipol, jpol, :,     ivar
                            utemp[ivar, jpol, ipol, :] = \
                                temp[:, np.argwhere(s_ids == ids[idx])[0][0]]
                ds_o[gll_idx] = utemp

    finally:
        try:
            f_in_1.close()
        except:
            pass
        try:
            f_in_2.close()
        except:
            pass
        try:
            f_in_3.close()
        except:
            pass
        try:
            f_in_4.close()
        except:
            pass
        try:
            f_out.close()
        except:
            pass


@click.command()
@click.argument("input_folder", type=click.Path(exists=True,
                                                file_okay=False,
                                                dir_okay=True))
@click.argument("output_folder", type=click.Path(exists=False))

def repack_database(input_folder, output_folder):
    found_filenames = []
    for root, _, filenames in os.walk(input_folder):
        for filename in filenames:
            if filename != "ordered_output.nc4":
                continue
            found_filenames.append(os.path.join(root, filename))

    assert found_filenames, "No files named `ordered_output.nc4` found."

    os.makedirs(output_folder)

    # The unrolled merge completely unrolls everything, dededuplicates the GLL
    # points, and merges both netCDF files into one big file.
    unroll_and_merge_netcdf4(filenames=found_filenames,
                             output_folder=output_folder)


if __name__ == "__main__":
    repack_database()
