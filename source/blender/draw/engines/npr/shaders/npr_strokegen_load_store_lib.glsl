#ifndef NPR_STROKEGEN_LOAD_STORE_LIB_INCLUDED
#define NPR_STROKEGEN_LOAD_STORE_LIB_INCLUDED


#ifdef CAT 
#   undef CAT
#endif

#ifdef CAT_
#   undef CAT_
#endif

/* Macro expansion, for details, see
/* ---------------------------------------
/* https://stackoverflow.com/questions/1489932/how-to-concatenate-twice-with-the-c-preprocessor-and-expand-a-macro-as-in-arg */
#define CAT(x, y) CAT_(x, y)
#define CAT_(x, y) x ## y


#define Load2(buf, index, val) \
    (val).x = buf[(index)*2+0]; \
    (val).y = buf[(index)*2+1]; \


#define Store2(buf, index, val) \
    buf[(index)*2+0] = (val).x; \
    buf[(index)*2+1] = (val).y; \


#define Load3(buf, index, val) \
    (val).x = buf[(index)*3+0]; \
    (val).y = buf[(index)*3+1]; \
    (val).z = buf[(index)*3+2]; \


#define Store3(buf, index, val) \
    buf[(index)*3+0] = (val).x; \
    buf[(index)*3+1] = (val).y; \
    buf[(index)*3+2] = (val).z; \


#define Load4(buf, index, val) \
    (val).x = buf[(index)*4+0]; \
    (val).y = buf[(index)*4+1]; \
    (val).z = buf[(index)*4+2]; \
    (val).w = buf[(index)*4+3]; \


#define Store4(buf, index, val) \
    buf[(index)*4+0] = (val).x; \
    buf[(index)*4+1] = (val).y; \
    buf[(index)*4+2] = (val).z; \
    buf[(index)*4+3] = (val).w; \







#endif


