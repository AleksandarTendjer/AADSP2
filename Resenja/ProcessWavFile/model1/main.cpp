/////////////////////////////////////////////////////////////////////////////////
//
// @file main.cpp
//
// Module: model0
// Description:  Add multitap echo to input signal
// $Source: $
// $Revision: 1.0 $
// $Date: 29.11.2018
// $Author:  Aleksandar Tendjer
//
/////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <string.h>
#include "stdfix_emu.h"
#include "WAVheader.h"
#include "tremolo1.h"
#include "common.h"


double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];
/////////////////////////////////////////////////////////////
//////GLOBALS
////////////////////////////////////////////////////////////////////////////////
tremolo_struct_t tremolo_data;
all_gains gains;
mode_sel selected=mode1; 
enable en=on;
double input_gain=(double)0.501187 ;
double headroom_gain=(double)0.707946;

double** globFirstCH ;
double** globSeccondCH ;
double** globThirdCH ;
double** globFourthCH ;
double** globFifthCH ;


/////////////////////////////////////////////////////////////////////////////////
// @Author	ALeksandar Tendjer
// @Date		29.11.2018
//
// Function:
// initialFun
//
// @param - echoState - Control state structure
//		  - buffer - buffer for keeping delayed samples
//		  - echoBufLen - Length of buffer
//		  - delay - array containing delay values in number of samples
//		  - input_gain - gain to be applayed to input sample
//		  - tap_gain - array of gains to be applayed to each delayed sample
//		  - n_tap - number of taps (equals length of delay and tap_gain)
//
// @return - nothing
// Comment: Initialize echoState structure
//
// E-mail:	<email>
//
/////////////////////////////////////////////////////////////////////////////////
void initFun( )
{
	gains.input_gain=input_gain;
	gains.headroom_gain=headroom_gain;
}
double dbToInt(double inNum)
{
	return pow(10, inNum / 10);
}
void proccessingFun()
{
	double* firstCH = &sampleBuffer[0][0];
	double* seccondCH = &sampleBuffer[1][0];
	double* thirdCH = &sampleBuffer[2][0];
	double* fourthCH = &sampleBuffer[3][0];
	double* fifthCH = &sampleBuffer[4][0];
	double* buffer;
	double temp;
	int i, k;

	for (k = 0; k < 2; k++)
	{
		//sam pokazivac povecavam,unutar prvog niza 
		buffer = &sampleBuffer[k][0];

		thirdCH = &sampleBuffer[2][0];
		for (i = 0; i < BLOCK_SIZE; i++)
		{
			//buffer = firstCH;
			*buffer = *buffer * gains.input_gain;
			temp = *buffer;
			*thirdCH += temp;
			//++firstCH;
			++buffer;
			++thirdCH;
		}	
	}
	if (selected == 2)
	{

		for (i = 0; i < BLOCK_SIZE; i++)
		{
			*(fourthCH+i) += *(firstCH+i);
			*(fifthCH+i) += *(seccondCH+i);
			/*++firstCH;
			++seccondCH;
			++fourthCH;
			++fifthCH;*/
		}
		processBlock(&sampleBuffer[3][0], &sampleBuffer[3][0], &tremolo_data, BLOCK_SIZE);
		processBlock(&sampleBuffer[4][0], &sampleBuffer[4][0], &tremolo_data, BLOCK_SIZE);
	}

	firstCH = &sampleBuffer[0][0];
	seccondCH = &sampleBuffer[1][0];
	thirdCH = &sampleBuffer[2][0];
		for (i = 0; i < BLOCK_SIZE; i++)
		{

			*thirdCH = *thirdCH*gains.headroom_gain;
			*firstCH = *thirdCH*gains.input_gain;
			*seccondCH = *thirdCH*gains.input_gain;
			++firstCH;
			++seccondCH;
			++thirdCH;
		}


}
//
//void proccessingFun(double** buffer)
//{
//	int i, k;
//
//	double** firstCH = buffer;
//	double** seccondCH = buffer + 1;
//	double** thirdCH = buffer + 2;
//	double** fourthCH = buffer + 3;
//	double** fifthCH = buffer + 4;
//	
//	double temp;
//	
//	for (k = 0; k < 2; k++)
//	{
//		for (i = 0; i < BLOCK_SIZE; i++)
//		{
//			*buffer = *firstCH ;
//			//*thirdCH+=i;
//			**buffer = **buffer * gains.input_gain;
//			temp = **buffer;		
//			*(*thirdCH) += temp;
//			++*firstCH;
//			++*thirdCH;
//		}
//		//sam pokazivac povecavam,unutar prvog niza  
//		buffer++;
//	}
//	
//	if (selected == 2)
//	{
//		buffer -= k;
//		 firstCH = buffer;
//		 seccondCH = buffer + 1;
//		 thirdCH = buffer + 2;
//		 fourthCH = buffer + 3;
//		 fifthCH = buffer + 4;
//
//		for (i = 0; i < BLOCK_SIZE; i++)
//		{
//			*(*fourthCH) += *(*firstCH);
//			*(*fifthCH) += *(*seccondCH);
//			++*firstCH;
//			++*seccondCH;
//			++*fourthCH;
//			++*fifthCH;
//		}
//		
//		processBlock(*fourthCH,*fourthCH,&tremolo_data,BLOCK_SIZE);
//		processBlock(*fifthCH, *fifthCH, &tremolo_data, BLOCK_SIZE);
//	}
//	thirdCH = buffer + 2;
//	seccondCH = buffer + 1;
//	for (i = 0; i < BLOCK_SIZE; i++)
//	{
//		
//		*(*thirdCH) = *(*thirdCH)*gains.headroom_gain;
//		*(*firstCH) = *(*thirdCH)*gains.input_gain;
//		*(*seccondCH) = *(*thirdCH)*gains.input_gain;
//		++*firstCH;
//		++*seccondCH;
//		++*thirdCH;
//	}
//}



























//void proccessingFun(double** buffer)
//{
//	double temp = 0.0;
//	//go through input buffers and do the first gain, it will done both the way
//	//dupli pokazivac vrednost 
//	double** readPtrCH ;
//	double** firstCH = buffer;
//	int i = 0;
//
//	//first two channels 
//	for (readPtrCH = firstCH; readPtrCH < firstCH + 2; readPtrCH++)
//	{
//		buffer = readPtrCH;
//		//16 members in the array
//		for (i = 0; i < 16;i++)//buffer=readPtrCH; *buffer < end; *buffer++)
//		{
//			*buffer = *readPtrCH;
//			**buffer =**readPtrCH*gains.input_gain;
//			temp = **buffer;
//
//			//pristuppamo trecem kanalu 
//			buffer = firstCH + 2;
//			**buffer += temp;			
//		}
//	}
//	//tremolo function
//	if (selected == 2)
//	{
//		for (readPtrCH = firstCH ; readPtrCH < firstCH +2; readPtrCH++)
//		{
//			//pristupamo 4 ili 5
//			temp = **readPtrCH;
//			buffer = readPtrCH + 3;
//			**buffer += temp;
//			processBlock(*buffer, *buffer, &tremolo_data, BLOCK_SIZE);
//		}
//	}
//
//	//doing a headroom gain and again the input gains
//		for(i=0;i<BLOCK_SIZE;i++)
//		{
//			buffer = firstCH;
//			**buffer = **readPtrCH*gains.input_gain;
//			//	buffer[0][i]=(buffer[2][i]* gains.input_gain);
//
//			buffer = firstCH + 1;
//			**buffer = **readPtrCH*gains.input_gain;
//			//buffer[1][i]=(buffer[2][i]* gains.input_gain);
//
//			buffer = firstCH + 2;
//			**buffer = **readPtrCH*gains.headroom_gain;
//				//buffer[2][i]=(buffer[2][i]* gains.headroom_gain);		
//		}			
//}
int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;
	
	// Init channel buffers
	for(int i=0; i<MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i],0,BLOCK_SIZE);
		
	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName,argv[1]);
	wav_in = OpenWavFileForRead (WavInputName,"rb");
	strcpy(WavOutputName,argv[2]);
	wav_out = OpenWavFileForRead (WavOutputName,"wb");
	//-------------------------------------------------
	// Read all arguments
	if(argv[3][0]!=NULL&&argv[3][0]=='0')
		en=off;
	/*if (argv[4][0] != NULL)
		input_gain = dbToInt(strtod(argv[4],NULL));
	if (argv[4][0] != NULL)
		headroom_gain = dbToInt(strtod(argv[4],NULL));*/
	if(argv[6][0] != NULL &&argv[6][0]=='1')
		selected=mode2;
	//-------------------------------------------------
	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	outputWAVhdr.fmt.NumChannels = MAX_NUM_CHANNEL;// inputWAVhdr.fmt.NumChannels; // change number of channels

	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size/inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate/inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign/inputWAVhdr.fmt.NumChannels;
	
	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out,outputWAVhdr);
	//initialize the tremolo data_structure and call the init function 
	init(&tremolo_data);
	initFun();
	
	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample / 8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size / (inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample / 8);

		//creating an array of all first channel element poiters

		double* buffPtr[] = {&sampleBuffer[0][0], &sampleBuffer[1][0], &sampleBuffer[2][0], &sampleBuffer[3][0], &sampleBuffer[4][0]};
		//printf(buffPtr);
		//pointer to an array of pointers 
		//double** buffPtrPtr = buffPtr;
		
		// exact file length should be handled correctly...
		for(int i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{	
			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<inputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);
					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;				// scale sample to 1.0/-1.0 range		
				}
			}
			
			
			//processing();
			proccessingFun();

			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<outputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = sampleBuffer[k][j] * SAMPLE_SCALE ;	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample/8, 1, wav_out);		
				}
			}		
		}
	}
	
	// Close files
	//-------------------------------------------------	
	fclose(wav_in);
	fclose(wav_out);
	//-------------------------------------------------	

	return 0;
}