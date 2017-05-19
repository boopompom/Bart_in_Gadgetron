/****************************************************************************************************************************
 * Description: Geometric Coil Compression Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Author: Jia Sen
 * Lang: C++
 * Date: 05/17/2017
 * Version: 0.0.1
****************************************************************************************************************************/

#ifndef BART_GCC_GADGET_H
#define BART_GCC_GADGET_H

#include "GenericReconGadget.h"
#include "GenericReconCartesianSpiritGadget.h"
#include "gadgetron_mricore_export.h"
#include "mri_core_data.h"
#include <gadgetron_paths.h>

#include "Bart_fileio.h"


#if defined (WIN32)
#ifdef __BUILD_GADGETRON_BartGccGadget__
#define EXPORTGADGETS_BartGccGadget __declspec(dllexport)
#else
#define EXPORTGADGETS_BartGccGadget __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_BartGccGadget
#endif


namespace Gadgetron {

	class EXPORTGADGETS_BartGccGadget BartGccGadget : public GenericReconDataBase
	{

	public:
		GADGET_DECLARE(BartGccGadget);
		
		typedef GenericReconDataBase BaseClass;

		BartGccGadget();
		virtual ~BartGccGadget() = default;

	protected:
		GADGET_PROPERTY(BartWorkingDirectory, std::string, "Absolute path to temporary file location (will default to workingDirectory)", "");
		GADGET_PROPERTY(BartWorkingDirectoryDelete, bool, "Whether to delete BartWorkingDirectory", true);
		
		GADGET_PROPERTY(CalibSize, int, "Size of CalibSize", 24);
		GADGET_PROPERTY(DstChaNum, int, "Compressed Channel Number",12);
                
		virtual int process_config(ACE_Message_Block* mb);
		virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
		
		long long image_counter_;
		std::string workLocation_;
		 
	};

}
#endif //BART_GCC_GADGET_H
