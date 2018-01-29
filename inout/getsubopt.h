#ifndef GETSUBOPT_H
#define GETSUBOPT_H

#ifdef __cplusplus
extern "C" {
#endif

    
int getsubopt_getsubopt(char **optionp,
			const char **tokens,
			char **valuep
			);

#ifdef __cplusplus
}
#endif 

#endif /*GETSUBOPT_H*/
