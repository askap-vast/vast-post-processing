from vast_post_processing import crop, corrections, compress

def run(data_root: Union[str, Path],
        crop_size: u.quantity.Quantity,
        epoch: Union[str, int, list],
        stokes: str,
        out_root: Optional[Union[str, Path]]=None,
        create_moc: Optional[bool]=False,
        overwrite: Optional[bool]=False,
        ):

    if out_root is None:
        out_root = data_root

    image_path_glob_list: list[Generator[Path, None, None]] = []
    
    image_root = data_root / f"STOKES{stokes}_IMAGES"
    logger.debug(image_root)
    
    if type(epoch) is int:
        epoch = list(epoch)
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(
            image_root.glob(f"epoch_*/*.fits")
        )
    else:
        for n in epoch:
            image_path_glob_list.append(
                image_root.glob(f"epoch_{n}/*.fits")
            )
    
    for image_path in chain.from_iterable(image_path_glob_list):
        logger.info(f"Working on {image_path}...")
        epoch_dir = image_path.parent.name
        _, _, field, sbid_str, *_ = image_path.name.split(".")
        sbid = int(sbid_str[2:])
        
        
        
        # get rms and background images
        rms_path = (
            data_root
            / f"STOKES{stokes}_RMSMAPS"
            / epoch_dir
            / f"noiseMap.{image_path.name}"
        )
        
        bkg_path = (
            data_root
            / f"STOKES{stokes}_RMSMAPS"
            / epoch_dir
            / f"meanMap.{image_path.name}"
        )
        
        # get selavy files
        components_name = f"selavy-{image_path.name}".replace(".fits",
                                                          ".components.xml"
                                                          )           
        islands_name = components_name.replace("components", "islands")
        
        selavy_dir = (
            data_root
            / f"STOKES{stokes}_SELAVY"
            / epoch_dir
        )
        components_path = selavy_dir / components_name
        islands_path = selavy_dir / islands_name
        
        exists = True
        if not rms_path.exists():
            exists = False
            logger.warning(f"noisemap file ({rms_path}) is missing.")
        
        if not bkg_path.exists():
            exists = False
            logger.warning(f"meanmap file ({bkg_path}) is missing.")
        if not components_path.exists():
            exists = False
            logger.warning(f"selavy components file ({components_path}) is missing.")
        if not islands_path.exists():
            exists = False
            logger.warning(f"selavy islands file ({islands_path}) is missing.")
        
        if not exists:
            # is best practice to skip or abort entirely?
            logger.warning(f"Skipping {image_path} due to missing files.")
        
        corrected_fits, corrected_cats = corrections.correct_field(image_path)
        
        # handle the fits files
        for i, path in enumerate((rms_path, bkg_path, image_path)):
            stokes_dir = f"{path.parent.parent.name}_CROPPED" # what suffix should we use?
            fits_output_dir = out_root / stokes_dir / epoch_dir
            if not fits_output_dir.exists():
                fits_output_dir.mkdir(parents=True)
            
            outfile = fits_output_dir / path.name
            hdu = corrected_fits[i]
            field_centre = get_field_centre(hdu.header)
            cropped_hdu = crop_hdu(hdu, field_centre, size=crop_size)
            cropped_hdu.writeto(outfile, overwrite=overwrite)
            logger.debug(f"Wrote {outfile}")
        
        
        # Crop the catalogues
        stokes_dir = f"{components_path.parent.parent.name}_CROPPED" # what suffix should we use?
        cat_output_dir = out_root / stokes_dir / epoch_dir
        
        if not cat_output_dir.exists():
            cat_output_dir.mkdir(parents=True)
        
        for i, path in enumerate((components_path, islands_path)):
            outfile = cat_output_dir / path.name
            vot = corrects_cats[i]
        
            # This uses the last cropped hdu from the previous for loop
            # which should be the image file, but doesn't actually matter
            cropped_vot = crop_catalogue(vot,
                                         cropped_hdu,
                                         field_centre,
                                         size
                                         )

            if outfile.exists() and not overwrite:
                logger.critical(f"{outfile} exists, not overwriting")
            else:
                vot.to_xml(str(outfile))
                logger.debug(f"Wrote {outfile}")
        
        # Create the MOC
        if create_moc:
            moc_dir = f"STOKES{stokes}_MOC_CROPPED"
            moc_output_dir = out_root / moc_dir / epoch_dir
            
            moc_filename = image_path.name.replace('.fits','.moc.fits')
            moc_outfile = moc_output_dir / moc_filename
            
            if not moc_output_dir.exists():
                moc_output_dir.mkdir(parents=True)
            moc = vpc.wcs_to_moc(cropped_hdu)
            moc.write(moc_outfile, overwrite=overwrite)
            logger.debug(f"Wrote {moc_outfile}")
            
            stmoc_filename = image_path.name.replace('.fits','.stmoc.fits')
            stmoc_outfile = moc_output_dir / stmoc_filename
            
            stmoc = vpc.moc_to_stmoc(moc, cropped_hdu)
            stmoc.write(stmoc_outfile, overwrite=overwrite)
            logger.debug("Wrote {stmoc_outfile}")
