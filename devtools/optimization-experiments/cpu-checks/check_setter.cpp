#define NOMINMAX
#include <windows.h>
#include <xmmintrin.h>
#include <openmm/Vec3.h>
#include <openmm/common/ComputeVectorTypes.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
using namespace std;
using namespace OpenMM;
static void require(bool ok, const char* why) { if (!ok) throw runtime_error(why); }
struct ThreadPool {
    struct Task { virtual void execute(ThreadPool&,int)=0; virtual ~Task() {} };
    int threads=4, mismatch=-1, calls=0, waits=0;
    int getNumThreads() const { return threads; }
    void execute(Task& task) {
        ++calls;
        // Execute the exact production partition/control logic in deterministic
        // mock workers. This does not claim to test native pool concurrency.
        unsigned control=_mm_getcsr();
        for(int i=0;i<threads;++i) {
            _mm_setcsr(i==mismatch ? control^0x2000 : control);
            task.execute(*this,i);
        }
        _mm_setcsr(control);
    }
    void waitForThreads() { ++waits; }
};
struct Array {
    int count, element, reportCount, reportElement, downloads=0, uploads=0;
    bool fail=false, tail=false;
    const void* prefix=nullptr;
    vector<unsigned char> bytes;
    vector<string>* events;
    string name;
    Array(int n,int e,vector<string>& ev,string label):count(n),element(e),reportCount(n),reportElement(e),bytes(n*e),events(&ev),name(label) {}
    int getSize() const { return reportCount; }
    int getElementSize() const { return reportElement; }
    void download(void* to) { ++downloads; events->push_back(name+".download"); memcpy(to,bytes.data(),bytes.size()); }
    void upload(const void* from) {
        ++uploads; events->push_back(name+".upload"); tail=from!=prefix;
        if(fail) throw runtime_error("injected upload failure");
        memcpy(bytes.data(),from,bytes.size());
    }
};
struct ComputeContext {
    int n,p,contexts=1,reorders=0;
    bool mixed=true,dbl=false;
    vector<string> events;
    Array posq,correction,velm;
    vector<int> order;
    vector<mm_int4> offsets;
    ThreadPool pool;
    void* allocation;
    unsigned char* pinned;
    ComputeContext(int count,bool useDouble=false,int misalign=0):n(count),p((count+31)/32*32),mixed(!useDouble),dbl(useDouble),
        posq(p,useDouble?32:16,events,"posq"),correction(p,16,events,"correction"),velm(p,32,events,"velm"),order(n),offsets(p) {
        allocation=_aligned_malloc(static_cast<size_t>(p)*32+64,32);
        require(allocation!=nullptr,"allocation"); pinned=(unsigned char*)allocation+misalign;
        memset(allocation,0xa5,static_cast<size_t>(p)*32+64);
        posq.prefix=correction.prefix=pinned;
        for(int i=0;i<n;++i) order[i]=n-i-1;
        for(int i=0;i<p;++i) {
            offsets[i]=mm_int4(1,2,3,4);
            if(dbl) {
                mm_double4 q(9,8,7,i%2?-0.0:0.375); memcpy(posq.bytes.data()+i*32,&q,32);
            } else {
                mm_float4 q(9,8,7,i%2?-0.0f:0.375f); memcpy(posq.bytes.data()+i*16,&q,16);
            }
        }
        memcpy(pinned,posq.bytes.data(),posq.bytes.size());
    }
    ~ComputeContext() { _aligned_free(allocation); }
    bool getUseMixedPrecision() const { return mixed; }
    bool getUseDoublePrecision() const { return dbl; }
    int getNumContexts() const { return contexts; }
    int getPaddedNumAtoms() const { return p; }
    Array& getPosq() { return posq; } Array& getPosqCorrection() { return correction; } Array& getVelm() { return velm; }
    ThreadPool& getThreadPool() { return pool; }
    void* getPinnedBuffer() { return pinned; }
    const vector<int>& getAtomIndex() { return order; }
    vector<mm_int4>& getPosCellOffsets() { return offsets; }
    void reorderAtoms() { ++reorders; events.push_back("reorder"); }
};
struct System { int n; int getNumParticles() { return n; } };
struct ContextImpl { System system; System& getSystem() { return system; } };
struct ContextSelector { explicit ContextSelector(ComputeContext&) {} };
struct Platform { string name="CUDA"; string getName() { return name; } };
// Execute the current upstream copyCoordinateBuffers body on CPU arrays.
// One mock worker models only its xyz write contract, not CUDA execution.
#define KERNEL
#define GLOBAL
#define RESTRICT
#define GLOBAL_ID 0
#define GLOBAL_SIZE 1
#define SUPPORTS_DOUBLE_PRECISION
using float4 = mm_float4;
using double4 = mm_double4;
#include "exact_copy_coordinates.inc"
#undef KERNEL
#undef GLOBAL
#undef RESTRICT
#undef GLOBAL_ID
#undef GLOBAL_SIZE
#undef SUPPORTS_DOUBLE_PRECISION
struct CopyKernel {
    Array& source;
    Array* dest=nullptr;
    bool useDouble, fail=false;
    CopyKernel(Array& source,bool useDouble):source(source),useDouble(useDouble) {}
    void setArg(int arg,Array& array) { require(arg==1,"copy argument"); dest=&array; }
    void execute(int count) {
        require(dest!=nullptr,"unbound copy");
        source.events->push_back(useDouble?"copyDoubleKernel.execute":"copyFloatKernel.execute");
        if(fail) throw runtime_error("injected copy failure");
        if(useDouble) {
            vector<double> in(source.count); vector<mm_double4> out(dest->count);
            memcpy(in.data(),source.bytes.data(),source.bytes.size());
            memcpy(out.data(),dest->bytes.data(),dest->bytes.size());
            copyDoubleBuffer(in.data(),out.data(),count);
            memcpy(dest->bytes.data(),out.data(),dest->bytes.size());
        } else {
            vector<float> in(source.count); vector<mm_float4> out(dest->count);
            memcpy(in.data(),source.bytes.data(),source.bytes.size());
            memcpy(out.data(),dest->bytes.data(),dest->bytes.size());
            copyFloatBuffer(in.data(),out.data(),count);
            memcpy(dest->bytes.data(),out.data(),dest->bytes.size());
        }
    }
};
struct CommonUpdateStateDataKernel {
    ComputeContext& cc; Platform platform;
    Array floatBuffer,doubleBuffer;
    CopyKernel floatKernel,doubleKernel;
    CopyKernel *copyFloatKernel,*copyDoubleKernel;
    explicit CommonUpdateStateDataKernel(ComputeContext& c):cc(c),
        floatBuffer(3*c.n,4,c.events,"floatBuffer"),doubleBuffer(3*c.n,8,c.events,"doubleBuffer"),
        floatKernel(floatBuffer,false),doubleKernel(doubleBuffer,true),
        copyFloatKernel(&floatKernel),copyDoubleKernel(&doubleKernel) {
        floatBuffer.prefix=doubleBuffer.prefix=c.pinned;
    }
    Platform& getPlatform() { return platform; }
    void setPositions(ContextImpl&,const vector<Vec3>&);
};
struct BaselineKernel : CommonUpdateStateDataKernel {
    using CommonUpdateStateDataKernel::CommonUpdateStateDataKernel;
    void setPositions(ContextImpl&,const vector<Vec3>&);
};
#include "exact_helper.inc"
#include "exact_candidate.inc"
#define CommonUpdateStateDataKernel BaselineKernel
#include "exact_baseline.inc"
#undef CommonUpdateStateDataKernel

static double bits(uint64_t value) { double d; memcpy(&d,&value,8); return d; }
static vector<Vec3> inputs(int n,bool safe) {
    vector<Vec3> v(n);
    vector<double> edges={0.0,-0.0,1.0,-1.0,1.000000059604644775390625,1.0000000596046449,
        numeric_limits<double>::denorm_min(),-numeric_limits<double>::denorm_min(),
        numeric_limits<float>::min()/2.0, numeric_limits<float>::denorm_min()/2.0,
        1e40,-1e40,numeric_limits<double>::max(),numeric_limits<double>::infinity(),
        -numeric_limits<double>::infinity(),bits(0x7ff8000000000123ULL),bits(0x7ff0000000000456ULL)};
    uint64_t seed=0x123456789abcULL;
    for(int i=0;i<n;++i) for(int j=0;j<3;++j) {
        seed=seed*6364136223846793005ULL+1442695040888963407ULL;
        v[i][j]=safe ? double((i+j)%4096) : ((i*3+j)%4 ? edges[(i*3+j)%edges.size()] : bits(seed));
    }
    return v;
}
static int cases=0, fusedCases=0, fallbackCases=0, failureCases=0;
static void run(int n,unsigned control,int variant=0,int mismatch=-1,int fail=0) {
    const unsigned saved=_mm_getcsr(); _mm_setcsr(control);
    bool dbl=variant==1;
    ComputeContext a(n,dbl,variant==7?4:0),b(n,dbl,variant==7?4:0);
    CommonUpdateStateDataKernel ca(a); BaselineKernel cb(b); ContextImpl context{{n}};
    a.pool.mismatch=b.pool.mismatch=mismatch;
    if(variant==2) a.mixed=b.mixed=false;
    if(variant==3) ca.platform.name=cb.platform.name="OpenCL";
    if(variant==4) a.contexts=b.contexts=2;
    if(variant==5) a.velm.reportElement=b.velm.reportElement=16;
    if(variant==6) a.posq.reportCount=b.posq.reportCount=a.p+1;
    if(variant==8) a.pool.threads=b.pool.threads=1;
    const char* flag=variant==10?"":variant==11?"0":variant==12?"10":variant==13?"true":"1";
    _putenv_s("OPENMM_EXPERIMENT_FUSED_SET_POSITIONS",flag);
    _putenv_s("OPENMM_EXPERIMENT_PARALLEL_SET_POSITIONS",variant==9?"0":"1");
    if(fail==1) ca.floatBuffer.fail=cb.floatBuffer.fail=ca.doubleBuffer.fail=cb.doubleBuffer.fail=true;
    if(fail==3) ca.floatKernel.fail=cb.floatKernel.fail=ca.doubleKernel.fail=cb.doubleKernel.fail=true;
    if(fail==2) a.correction.fail=b.correction.fail=true;
    bool safe=(control&0x1f80)!=0x1f80;
    vector<Vec3> v=inputs(n,safe),before=v;
    const auto charges=a.posq.bytes;
    string ae,be;
    try { ca.setPositions(context,v); } catch(const runtime_error& e) { ae=e.what(); }
    _mm_setcsr(control);
    try { cb.setPositions(context,v); } catch(const runtime_error& e) { be=e.what(); }
    _mm_setcsr(saved);
    require(ae==be,"exception difference");
    require(a.posq.bytes==b.posq.bytes,"base float bytes differ");
    require(a.correction.bytes==b.correction.bytes,"correction float bytes differ");
    require(a.events==b.events,"download/upload/reorder sequence differs");
    require(a.reorders==b.reorders,"reorder count differs");
    require(memcmp(a.offsets.data(),b.offsets.data(),a.offsets.size()*sizeof(mm_int4))==0,"offsets differ");
    require(memcmp(v.data(),before.data(),v.size()*sizeof(Vec3))==0,"input modified");
    for(int i=0;i<n;++i) require(memcmp(a.posq.bytes.data()+i*a.posq.element+3*a.posq.element/4,charges.data()+i*a.posq.element+3*a.posq.element/4,a.posq.element/4)==0,"charge changed");
    const bool fused=n>=4096 && !safe && (variant==0||variant==8||variant==9);
    if(fused && (fail==0 || fail==2)) require(a.correction.tail,"eligible branch did not use tail");
    if(!fused && a.mixed && fail!=1) require(!a.correction.tail,"ineligible branch used tail");
    require(a.posq.downloads==0 && b.posq.downloads==0,"packed setter must not download posq");
    for(int i=n;i<a.p;++i) require(memcmp(a.posq.bytes.data()+i*a.posq.element,charges.data()+i*a.posq.element,a.posq.element)==0,"padding posq changed");
    if(fail) { require(a.reorders==0,"reorder after failed upload"); ++failureCases; }
    else require(a.reorders==1,"missing reorder");
    if(fused) ++fusedCases; else ++fallbackCases;
    ++cases;
}
int main() {
    try {
        for(unsigned rounding : {0u,0x2000u,0x4000u,0x6000u}) for(unsigned denorm : {0u,0x40u,0x8000u,0x8040u})
            for(int n : {0,1,4095,4096,4097,8193}) for(int mismatch : {-1,0,2})
                run(n,0x1f80|rounding|denorm,0,mismatch);
        for(int variant=1;variant<=13;++variant) run(4097,0x1f80,variant,2);
        for(unsigned mask : {0x1f00u,0x1e80u,0x1d80u,0x1b80u,0x1780u,0x0f80u}) run(4097,mask);
        for(int fail : {1,2,3}) for(int mismatch : {-1,2}) run(4097,0x1f80,0,mismatch,fail);
        run(919221,0x1f80,0,2);
        cout << "{\"cases\":"<<cases<<",\"fused_cases\":"<<fusedCases<<",\"fallback_cases\":"<<fallbackCases<<",\"failure_cases\":"<<failureCases<<",\"exact_output_bytes\":true,\"mock_workers_only\":true,\"GPU_executed\":false}" << endl;
        return 0;
    } catch(const exception& e) { cerr<<e.what()<<endl; return 1; }
}
