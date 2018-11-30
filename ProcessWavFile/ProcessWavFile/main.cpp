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
#include "WAVheader.h"
#include "tremolo1.h"


/////////////////////////////////////////////////////////////////////////////////
///CONSTANTS
/////////////////////////////////////////////////////////////////////////////////
#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 5
#define N_TAP 4
//input gain is -6 db  is 0.5 
//power ratio is 0.25
//input gain is --3 db  is 0.25 
//power ratio is 0.25

double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];

typedef struct  
{
	double input_gain;
	double headroom_gain;
}all_gains;
typedef enum 
{
	mode1=1,
	mode2=2
}mode_sel;


/////////////////////////////////////////////////////////////
//////GLOBALS
////////////////////////////////////////////////////////////////////////////////
tremolo_struct_t tremolo_data;
all_gains al_gains;
mode_sel selected; 
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
void initFun(all_gains* gains,double headroom_gain,double input_gain)
{
	
		gains->input_gain=input_gain;
		gains->headroom_gain=headroom_gain;
}

void proccessingFun(double buffer[][BLOCK_SIZE],all_gains* gains)
{

	int i,k;
	//go through input buffers and do the first gain, it will done both the way
	for(k=0;k<2;k++)
		for(i=0;i<BLOCK_SIZE;i++)
		{
			buffer[k][i]=buffer[k][i]*gains->input_gain;
				buffer[2][i]+=buffer[k][i];
		}
	//tremolo function
		if(selected==2)
			for(k=3;k<5;k++)
			{
				processBlock(buffer[k],buffer[k],&tremolo_data,BLOCK_SIZE);
			}
	//doing a headroom gain and again the input gains
		for(i=0;i<BLOCK_SIZE;i++)
		{
			buffer[2][i]=(buffer[2][i]*gains->headroom_gain);
			//dal += ili =
			buffer[0][i]=(buffer[2][i]*gains->input_gain);
			buffer[1][i]=(buffer[2][i]*gains->input_gain);
			if(selected==2)
				buffer[3][i]+=(buffer[2][i]*gains->input_gain);
		}			
}
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

	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	outputWAVhdr.fmt.NumChannels = inputWAVhdr.fmt.NumChannels; // change number of channels

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
	initFun(&al_gains,0.707946,0.501187);
	
	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample/8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size/(inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample/8);
		
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
			selected=mode1;
			proccessingFun(sampleBuffer,&al_gains);

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