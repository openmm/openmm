/*
   This file is an interface between TINKER Fortran code
   and a Java Server used for socket based communication.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jni.h>

/* Rename functions for Intel Compilers on Windows */

#ifdef windowsintel
#define createjvm_ CREATEJVM
#define destroyjvm_ DESTROYJVM
#define getmonitor_ GETMONITOR
#define releasemonitor_ RELEASEMONITOR
#define createserver_ CREATESERVER
#define destroyserver_ DESTROYSERVER
#define needupdate_ NEEDUPDATE
#define setupdated_ SETUPDATED
#define createsystem_ CREATESYSTEM
#define createupdate_ CREATEUPDATE
#define setcoordinates_ SETCOORDINATES
#define setconnectivity_ SETCONNECTIVITY
#define setvdw_ SETVDW
#define setcharge_ SETCHARGE
#define setmass_ SETMASS
#define setatomic_ SETATOMIC
#define setatomtypes_ SETATOMTYPES
#define setname_ SETNAME
#define setstory_ SETSTORY
#define setkeyword_ SETKEYWORD
#define setforcefield_ SETFORCEFIELD
#define settime_ SETTIME
#define setenergy_ SETENERGY
#define setstep_ SETSTEP
#define settype_ SETTYPE
#define setfile_ SETFILE
#define setvelocity_ SETVELOCITY
#define setacceleration_ SETACCELERATION
#define setinduced_ SETINDUCED
#define setgradients_ SETGRADIENTS
#endif

/* Some global variables */

JNIEnv* env = 0;
jobject serverobject = 0;
jclass  serverclass = 0;
jobject systemobject = 0;
jclass  systemclass = 0;
jobject updateobject = 0;
jclass  updateclass = 0;

jboolean InitializeJVM() {
   JavaVMInitArgs args;
   JavaVMOption options[1];
   JavaVM *vm;
   jint numOptions = 1;
   jint r;
   char *classpath;
   char *def;

   classpath = getenv("CLASSPATH");
   if (classpath == NULL) {
      return JNI_FALSE;
   }
   //printf("%s\n", classpath);

   def = (char *) malloc(strlen(classpath) + 20);
   sprintf(def, "-Djava.class.path=%s", classpath);
   options[0].optionString = def;

   args.version = JNI_VERSION_1_6;
   args.nOptions = numOptions;
   args.options = options;
   args.ignoreUnrecognized = JNI_TRUE;

   r = JNI_CreateJavaVM(&vm, (void **) &env, &args);
   if (r == JNI_OK) {
      //printf("JNI_CreateJavaVM Success: %i\n", r);
      return JNI_TRUE;
   } else {
      //printf("JNI_CreateJavaVM Error: %i\n", r);
      return JNI_FALSE;
   }
}

void chksocket_(int *flag) {
   *flag = 1;
   return;
}

void createjvm_(int *flag) {
   if (!InitializeJVM()) {
      *flag = 0;
      return;
   }
   *flag = 1;
   return;
}

void destroyjvm_() {
   JavaVM *vm;
   (*env)->GetJavaVM(env, &vm);
   (*vm)->DestroyJavaVM(vm);
}

void getmonitor_() {
   (*env)->MonitorEnter(env, updateobject);
}

void releasemonitor_() {
   (*env)->MonitorExit(env, updateobject);
}

void createserver_(jint *flag) {
   jmethodID methodID;
   jobject temp;

   // Find the Java ffe.tinker.TinkerServer class
   serverclass = (*env)->FindClass(env, "ffe/tinker/TinkerServer");
   if (serverclass == 0) {
      printf("Could not find ffe.tinker.TinkerServer\n");
      *flag = 0;
      return;
   }
   methodID = (*env)->GetMethodID(env, serverclass, "<init>",
      "(Lffe/tinker/TinkerSystem;)V");
   temp =  (*env)->NewObject(env, serverclass, methodID, systemobject);
   methodID = (*env)->GetMethodID(env, serverclass, "start", "()V");
   (*env)->CallVoidMethod(env, temp, methodID, "()V");
   serverobject =  (*env)->NewGlobalRef(env, temp);
   (*env)->DeleteLocalRef(env, temp);
   *flag = 1;
   return;
}

void destroyserver_() {
   jmethodID methodID;
   jboolean jffe = JNI_FALSE;

   // Tell the server to close down
   methodID = (*env)->GetMethodID(env, serverclass, "stop", "()V");
   (*env)->CallVoidMethod(env, serverobject, methodID, "()V");

   // Wait while it closes any clients
   methodID = (*env)->GetMethodID(env, serverclass, "isAlive", "()Z");
   while ((*env)->CallBooleanMethod(env, serverobject, methodID, "()Z")) {
   }

   (*env)->DeleteGlobalRef(env, serverobject);
   return;
}

void needupdate_(int *status) {
   jmethodID methodID;
   jboolean ret;

   methodID = (*env)->GetMethodID(env, serverclass, "needUpdate", "()Z");
   ret = (*env)->CallBooleanMethod(env, serverobject, methodID, "()Z");
   if (ret) *status = 1;
   else *status = 0;
   return;
}

void setupdated_() {
   jmethodID methodID;

   methodID = (*env)->GetMethodID(env, serverclass, "setUpdated", "()V");
   (*env)->CallVoidMethod(env, serverobject, methodID, "()V");
   return;
}

jstring char2jstring(const char *str, int len) {
   jstring result;
   jbyteArray bytes = 0;
   jmethodID methodID;
   jclass class;

   if ((*env)->EnsureLocalCapacity(env, 2) < 0) {
      printf("Out of memory\n");
      exit(-1);
   }
   bytes = (*env)->NewByteArray(env, len);
   if (bytes != NULL) {
      (*env)->SetByteArrayRegion(env, bytes, 0, len, (jbyte *)str);
      class = (*env)->FindClass(env, "Ljava/lang/String;");
      methodID = (*env)->GetMethodID(env, class, "<init>", "([B)V");
      result = (*env)->NewObject(env, class, methodID, bytes);
      (*env)->DeleteLocalRef(env, bytes);
      if (result != NULL) return result;
    } /* else fall through */
    printf("Problem creating java string for: %s\nLength: %i", str, len);
    exit(-1);
}

void createsystem_(jint *atoms, jint *keywords, jint* flag) {
   jmethodID methodID;
   jobject temp;

   systemclass = (*env)->FindClass(env, "ffe/tinker/TinkerSystem");
   if (systemclass == 0) {
      printf("Could not find ffe.tinker.TinkerSystem\n");
      *flag = 0;
      return;
   }
   methodID = (*env)->GetMethodID(env, systemclass, "<init>", "(II)V");
   temp = (*env)->NewObject(env, systemclass, methodID, *atoms, *keywords);
   systemobject = (*env)->NewGlobalRef(env, temp);
   (*env)->DeleteLocalRef(env, temp);
   *flag = 1;
   return;
}

void createupdate_(jint *n, jint *type, jint* amoeba, jint* flag) {
   jmethodID methodID;
   jobject temp;
   jboolean jbool;

   updateclass = (*env)->FindClass(env, "ffe/tinker/TinkerUpdate");
   if (updateclass == 0) {
      printf("Could not find ffe.tinker.TinkerUpdate\n");
      *flag = 0;
      return;
   }

   methodID = (*env)->GetMethodID(env, updateclass, "<init>", "(IIZ)V");

   if (amoeba == 0) jbool = JNI_FALSE;
   else jbool = JNI_TRUE;

   temp = (*env)->NewObject(env, updateclass, methodID, *n, *type, jbool);
   updateobject  = (*env)->NewGlobalRef(env, temp);
   (*env)->DeleteLocalRef(env, temp);
   methodID = (*env)->GetMethodID(env,serverclass,"loadUpdate",
                                   "(Lffe/tinker/TinkerUpdate;)V");
   (*env)->CallObjectMethod(env, serverobject, methodID, updateobject);
   *flag = 1;
   return;
}

void setcoordinates_(jint *n, jdouble *x, jdouble *y, jdouble *z) {
   jfieldID fieldID;
   jobjectArray coords;
   jdoubleArray jx;
   jobject object;
   jclass class;

   object = updateobject;
   class = updateclass;
   if (updateobject == 0) {
      object = systemobject;
      class = systemclass;
   }
   fieldID = (*env)->GetFieldID(env, class, "coordinates", "[[D");
   coords = (*env)->GetObjectField(env, object, fieldID);
   jx = (*env)->GetObjectArrayElement(env, coords, 0);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, x);
   jx = (*env)->GetObjectArrayElement(env, coords, 1);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, y);
   jx = (*env)->GetObjectArrayElement(env, coords, 2);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, z);
}

void setconnectivity_(jint *n, jint *a, jint *b, jint *c, jint *d) {
   jfieldID fieldID;
   jobjectArray connect;
   jintArray ci;

   fieldID = (*env)->GetFieldID(env, systemclass, "connectivity", "[[I");
   connect = (*env)->GetObjectField(env, systemobject, fieldID);
   ci = (*env)->GetObjectArrayElement(env, connect, 0);
   (*env)->SetIntArrayRegion(env, ci, 0, *n, a);
   ci = (*env)->GetObjectArrayElement(env, connect, 1);
   (*env)->SetIntArrayRegion(env, ci, 0, *n, b);
   ci = (*env)->GetObjectArrayElement(env, connect, 2);
   (*env)->SetIntArrayRegion(env, ci, 0, *n, c);
   ci = (*env)->GetObjectArrayElement(env, connect, 3);
   (*env)->SetIntArrayRegion(env, ci, 0, *n, d);
}

/*
void setvdw_(jint *n, jdouble *vdw) {
   jfieldID fieldID;
   jdoubleArray jvdw;

   fieldID = (*env)->GetFieldID(env, systemclass, "vdw", "[D");
   jvdw = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetDoubleArrayRegion(env, jvdw, 0, *n, *vdw);
}
*/

void setcharge_(jint *n, jdouble *charge) {
   jfieldID fieldID;
   jdoubleArray jcharge;

   fieldID = (*env)->GetFieldID(env, systemclass, "charge", "[D");
   jcharge = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetDoubleArrayRegion(env, jcharge, 0, *n, charge);
}

void setmass_(jint *n, jdouble *mass) {
   jfieldID fieldID;
   jdoubleArray jmass;

   fieldID = (*env)->GetFieldID(env, systemclass, "mass", "[D");
   jmass = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetDoubleArrayRegion(env, jmass, 0, *n, mass);
}

void setatomic_(jint *n, jint *atomic) {
   jfieldID fieldID;
   jintArray jatomic;

   fieldID = (*env)->GetFieldID(env, systemclass, "atomic", "[I");
   jatomic = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetIntArrayRegion(env, jatomic, 0, *n, atomic);
}

void setatomtypes_(jint *n, jint *t) {
   jfieldID fieldID;
   jintArray types;

   fieldID = (*env)->GetFieldID(env, systemclass, "types", "[I");
   types = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetIntArrayRegion(env, types, 0, *n, t);
}

void setname_(jint *num, char *s, jint len) {
   jfieldID fieldID;
   jstring string;
   jobjectArray name;
   jint index;

   index = (*num) - 1;
   if (index < 0) {
      printf("Negative index into name array\n");
      exit(-1);
   }
    string = char2jstring(s, len);
   fieldID = (*env)->GetFieldID(env, systemclass, "name",
                                "[Ljava/lang/String;");
   name = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetObjectArrayElement(env, name, index, string);
}

void setstory_(jint *num, char *s, jint len) {
   jfieldID fieldID;
   jstring string;
   jobjectArray story;
   jint index;

   index = (*num) - 1;
   if (index < 0) {
      printf("Negative index into story array\n");
      exit(-1);
   }
   string = char2jstring(s, len);
   fieldID = (*env)->GetFieldID(env, systemclass, "story",
                                "[Ljava/lang/String;");
   story = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetObjectArrayElement(env, story, index, string);
}

void setkeyword_(jint *num, char *keyword, jint len) {
   jfieldID fieldID;
   jstring string;
   jobjectArray keywords;
   jint index;

   index = (*num) - 1;
   if (index < 0) {
      printf("Negative index into keyword array\n");
      exit(0);
   }

   string = char2jstring(keyword, len);
   fieldID = (*env)->GetFieldID(env, systemclass, "keywords",
                                "[Ljava/lang/String;");
   keywords = (*env)->GetObjectField(env, systemobject, fieldID);
   (*env)->SetObjectArrayElement(env, keywords, index, string);
}

void setforcefield_(char *forcefield, jint len) {
   jfieldID fieldID;
   jstring string;

   string = char2jstring(forcefield, len);
   fieldID = (*env)->GetFieldID(env, systemclass, "forcefield",
                                "Ljava/lang/String;");
   (*env)->SetObjectField(env, systemobject, fieldID, string);
}

void setmdtime_(jdouble *time) {
   jfieldID fieldID;

   fieldID = (*env)->GetFieldID(env, updateclass, "time", "D");
   (*env)->SetDoubleField(env, updateobject, fieldID, *time);
}

void setenergy_(jdouble *e) {
   jfieldID fieldID;

   fieldID = (*env)->GetFieldID(env, updateclass, "energy", "D");
   (*env)->SetDoubleField(env, updateobject, fieldID, *e);
}

void setstep_(jint *step) {
   jfieldID fieldID;

   fieldID = (*env)->GetFieldID(env, updateclass, "step", "I");
   (*env)->SetIntField(env, updateobject, fieldID, *step);
}

void settype_(jint *type) {
   jfieldID fieldID;

   fieldID = (*env)->GetFieldID(env, updateclass, "type", "I");
   (*env)->SetIntField(env, updateobject, fieldID, *type);
}

void setfile_(char *name, jint len) {
   jfieldID fieldID;
   jstring string;

   string = char2jstring(name, len);
   fieldID = (*env)->GetFieldID(env, systemclass, "file",
                                "Ljava/lang/String;");
   (*env)->SetObjectField(env, systemobject, fieldID, string);
}

void setvelocity_(jint *n, jdouble *x, jdouble *y, jdouble *z) {
   jfieldID fieldID;
   jobjectArray ind;
   jdoubleArray jx;

   fieldID = (*env)->GetFieldID(env, updateclass, "velocity", "[[D");
   ind = (*env)->GetObjectField(env, updateobject, fieldID);
   jx = (*env)->GetObjectArrayElement(env, ind, 0);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, x);
   jx = (*env)->GetObjectArrayElement(env, ind, 1);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, y);
   jx = (*env)->GetObjectArrayElement(env, ind, 2);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, z);
}

void setacceleration_(jint *n, jdouble *x, jdouble *y, jdouble *z) {
   jfieldID fieldID;
   jobjectArray ind;
   jdoubleArray jx;

   fieldID = (*env)->GetFieldID(env, updateclass, "acceleration", "[[D");
   ind = (*env)->GetObjectField(env, updateobject, fieldID);
   jx = (*env)->GetObjectArrayElement(env, ind, 0);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, x);
   jx = (*env)->GetObjectArrayElement(env, ind, 1);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, y);
   jx = (*env)->GetObjectArrayElement(env, ind, 2);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, z);
}

void setinduced_(jint *n, jdouble *x, jdouble *y, jdouble *z) {
   jfieldID fieldID;
   jobjectArray ind;
   jdoubleArray jx;

   fieldID = (*env)->GetFieldID(env, updateclass, "induced", "[[D");
   ind = (*env)->GetObjectField(env, updateobject, fieldID);
   jx = (*env)->GetObjectArrayElement(env, ind, 0);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, x);
   jx = (*env)->GetObjectArrayElement(env, ind, 1);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, y);
   jx = (*env)->GetObjectArrayElement(env, ind, 2);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, z);
}

void setgradients_(jint *n, jdouble *x, jdouble *y, jdouble *z) {
   jfieldID fieldID;
   jobjectArray ind;
   jdoubleArray jx;

   fieldID = (*env)->GetFieldID(env, updateclass, "gradients", "[[D");
   ind = (*env)->GetObjectField(env, updateobject, fieldID);
   jx = (*env)->GetObjectArrayElement(env, ind, 0);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, x);
   jx = (*env)->GetObjectArrayElement(env, ind, 1);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, y);
   jx = (*env)->GetObjectArrayElement(env, ind, 2);
   (*env)->SetDoubleArrayRegion(env, jx, 0, *n, z);
}
