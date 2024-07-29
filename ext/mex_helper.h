#include "mex.hpp"
#include "mexAdapter.hpp"

using matlab::mex::ArgumentList;
using matlab::engine::MATLABEngine;
using matlab::data::ArrayFactory;
using matlab::data::Array;

class MexHelper : public matlab::mex::Function {
protected:
    ArrayFactory factory;
    std::shared_ptr<MATLABEngine> matlabPtr;
public:
    MexHelper() {
        matlabPtr = getEngine();
    }

    // Helper functions to display messages in the MATLAB command prompt.
    void displayMessage(std::string message) {
        matlabPtr->feval(u"fprintf", 0, std::vector<Array>({ factory.createScalar(message) }));
    }

    void displayError(std::string errorMessage)
    {
        matlabPtr->feval(u"error", 0, std::vector<Array>({ factory.createScalar(errorMessage) }));
    }

    class mexstream : public std::streambuf {
    private:
        MexHelper&                        func;
        ArrayFactory&                    factory;
        std::shared_ptr<MATLABEngine>  matlabPtr;
    public:
        mexstream(MexHelper& _func,
                  std::shared_ptr<MATLABEngine>& _matlabPtr,
                  ArrayFactory& _factory)
        : func(_func), matlabPtr(_matlabPtr), factory(_factory) {}
    protected:
        virtual std::streamsize xsputn(const char *s, std::streamsize n) { func.displayMessage(std::string(s, n)); return n; }
        virtual int overflow(int c=EOF) { if (c != EOF) { func.displayMessage(std::string(1, static_cast<char>(c))); } return c; }
    };
    
    class scoped_redirect_cout {
    public:
        scoped_redirect_cout(MexHelper& _func,
                             std::shared_ptr<MATLABEngine>& _matlabPtr,
                             ArrayFactory& _factory)
        : mout(_func, _matlabPtr, _factory) { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
        ~scoped_redirect_cout() { std::cout.rdbuf(old_buf); }
    private:
        mexstream mout;
        std::streambuf *old_buf;
    };
};