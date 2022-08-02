from plantcv import plantcv as pcv
import os
import sys, traceback
import cv2
import numpy as np
import argparse
import string

def options():
    parser = argparse.ArgumentParser(description="Imaging processing with opencv")
    parser.add_argument("-i", "--image", help="Input image file.", required=True)
    parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=False)
    parser.add_argument("-r","--result", help="result file.", required= False )
    parser.add_argument("-r2","--coresult", help="result file.", required= False )
    parser.add_argument("-w","--writeimg", help="write out images.", default=False)
    parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", default=None)
    args = parser.parse_args()
    return args
     
def main():
    args = options()
    pcv.params.debug = args.debug
    pcv.params.debug_outdir = args.outdir

    img, path, filename = pcv.readimage(filename=args.image)

    img_wb = pcv.white_balance(img, mode='max', roi=(600, 250, 80, 80))

    s = pcv.rgb2gray_hsv(rgb_img=img_wb, channel='s')

    s_thresh = pcv.threshold.binary(gray_img=s, threshold=30, max_value=255, object_type='light')

    b = pcv.rgb2gray_lab(rgb_img=img_wb, channel='b')

    b_thresh = pcv.threshold.binary(gray_img=b, threshold=133, max_value=255, object_type='light')

    bs = pcv.logical_and(bin_img1=s_thresh, bin_img2=b_thresh)

    fill = pcv.fill(bs, 50)

    masked = pcv.apply_mask(img=img_wb, mask=fill, mask_color='white')

    id_objects, obj_hierarchy = pcv.find_objects(img=masked, mask=fill)

    roi_contour, roi_hierarchy= pcv.roi.rectangle(img=img_wb, x=850, y=1100, h=620, w=700)

    roi_objects, hierarchy, kept_mask, obj_area = pcv.roi_objects(img_wb, roi_contour, roi_hierarchy, id_objects, obj_hierarchy, 'partial')

    pcv.params.line_thickness = 5

    obj, mask = pcv.object_composition(img=img_wb, contours=roi_objects, hierarchy=hierarchy)

    shape_img = pcv.analyze_object(img=img_wb, obj=obj, mask=mask)

    color_histogram = pcv.analyze_color(img_wb, mask, 'hsv')
    
    hue_circular_mean = pcv.outputs.observations['hue_circular_mean']['value']

    height_img = pcv.analyze_bound_horizontal(img=img_wb, obj=obj, mask=mask, line_position=1755)

    pseudocolored_img = pcv.visualize.pseudocolor(gray_img=s, mask=kept_mask, cmap='jet')

    a = pcv.rgb2gray_hsv(img_wb,"v")
    inv_mask = pcv.invert(mask)
    pmasked = cv2.applyColorMap(a,cv2.COLORMAP_JET)
    mask1 = pcv.apply_mask(mask=mask,img=pmasked,mask_color="black")
    masked1 = pcv.apply_mask(mask=inv_mask,img=img_wb,mask_color="black")
    sum_img = pcv.image_add(mask1, masked1)

    iy,ix,iz = np.shape(img_wb)
    blank = (np.ones((iy,500,3),dtype=np.uint8))*255
    vis1 = np.concatenate((img_wb, blank), axis=1)
    vis2 = np.concatenate((vis1, sum_img), axis=1)
    vis3 = cv2.resize(vis2, (720, 251))
    outname = args.outdir+"/"+str(filename[:-4])+"_pseudocolor.png"
    pcv.print_image(vis3,outname)
    
    masked1, binary, contours, hierarchy = pcv.rectangle_mask(img=masked, p1=(0,0), p2=(550,2055), color="black")
    
    masked2, binary, contours, hierarchy = pcv.rectangle_mask(img=masked1, p1=(2000,0), p2=(2500,2055), color="black")
    
    nbmask = pcv.naive_bayes_classifier(masked2, "brachy_naive_bayes_pdfs.txt")
    
    classified_img = pcv.visualize.colorize_masks(masks=[nbmask['Healthy'], nbmask['Unhealthy'], nbmask['Background']], colors=['lime green', 'red', 'gray'])
    
    unhealthy_plant = np.count_nonzero(nbmask['Unhealthy'])
    healthy_plant = np.count_nonzero(nbmask['Healthy'])
    percent_unhealthy = unhealthy_plant / (unhealthy_plant + healthy_plant)

    pcv.outputs.add_observation(variable='percent_unhealthy', trait='percent of plant detected to be unhealthy', method='ratio of pixels', scale='percent', datatype=float,value=percent_unhealthy, label='percent')

    pcv.print_results(filename=args.result)

    pcv.outputs.clear()
    
    if args.coresult is not None:
        nirpath = pcv.get_nir(path,filename)
        nir, path1, filename1 = pcv.readimage(nirpath)
        nir2 = cv2.imread(nirpath,0)

        nmask = pcv.resize(img=mask, resize_x=0.13, resize_y=0.125)

        newmask = pcv.crop_position_mask(img=nir, mask=nmask, x=20, y=3, v_pos="top", h_pos="left")

        nir_objects, nir_hierarchy = pcv.find_objects(img=nir, mask=newmask)

        pcv.params.line_thickness = 1

        nir_combined, nir_combinedmask = pcv.object_composition(img=nir, contours=nir_objects, hierarchy=nir_hierarchy)

        nir_hist = pcv.analyze_nir_intensity(gray_img=nir2, mask=nir_combinedmask, bins=256, histplot=True)

        nir_shape_image = pcv.analyze_object(img=nir2, obj=nir_combined, mask=nir_combinedmask)

        pcv.print_image(img=nir_hist, filename=os.path.join(pcv.params.debug_outdir, 'nirhist.png'))
        pcv.print_image(img=nir_shape_image, filename=os.path.join(pcv.params.debug_outdir, 'nirshape.png'))

        pcv.print_results(filename=args.coresult)
        
    pcv.outputs.clear()

if __name__ == '__main__':
    main()

