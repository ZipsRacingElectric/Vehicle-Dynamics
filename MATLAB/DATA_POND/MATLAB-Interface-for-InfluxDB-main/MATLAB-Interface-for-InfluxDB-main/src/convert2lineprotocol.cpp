#include <mex.hpp>
#include <mexAdapter.hpp>
#include "MatlabDataArray/ArrayType.hpp"
#include <vector>
#include <string>
#include <cmath>

//# Copyright 2024 - 2025 The MathWorks, Inc.

using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

        TypedArray<MATLABString> measurement = std::move(inputs[0]);
        TypedArray<MATLABString> fieldNames = std::move(inputs[1]);
        TypedArray<MATLABString> tagNames = std::move(inputs[2]);
        CellArray data = std::move(inputs[3]);
        uint64_t startRow = inputs[4][0];
        uint64_t endRow = inputs[5][0];

        std::u16string res = u"";

        std::u16string measurementStr = *(measurement.begin());

        TypedArrayRef<int64_t> timeDataArray = data[0];
        
        bool isTagNamesEmpty = false;
        if(tagNames.end()-tagNames.begin()<=0) isTagNamesEmpty = true;

        auto lastTag = tagNames.end();
        --lastTag;

        auto lastMinusOne = fieldNames.end() - 1;

        for (uint64_t iRow = startRow; iRow < endRow; iRow++) {
            std::u16string rowStr = u"";
            // Flag to indicate if all the field in the row is NaNs. If they are, addRow = false.
            bool addRow = true;
            int nanNum = 0;
            rowStr += measurementStr;
            if(isTagNamesEmpty) rowStr += u" ";
            else rowStr += u",";
               
            int colIndx = 0;

            for (auto it = tagNames.begin(); it != tagNames.end(); ++it) {
                rowStr += static_cast<std::u16string>(*it) + u"=";
                colIndx++;
                ArrayRef col = data[colIndx];
                switch (col.getType()) {
                    case ArrayType::CELL: {
                        // cell array of char arrays. g3432359
                        CharArrayRef entry = col[iRow][0];
                        String strEntry = entry.toUTF16();
                        rowStr += strEntry;
                        break;
                    }
                    case ArrayType::MATLAB_STRING: {
                         TypedArrayRef<MATLABString> col = data[colIndx];
                        std::u16string entry = col[iRow];
                        rowStr += entry;
                        break;
                    }
                    default: {
                        break;
                    }

                }
               
                if (it == lastTag) rowStr += u" ";
                else rowStr += u",";
            }

            for (auto it = fieldNames.begin(); it != fieldNames.end(); ++it) {
                std::u16string temp = static_cast<std::u16string>(*it) + u"=";
                colIndx++;
                bool is_nan = false;
                ArrayRef col = data[colIndx];

                switch (col.getType())
                {
                    case ArrayType::INT8: {
                        int8_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"i";
                        break;
                    }
                    case ArrayType::UINT8: {
                        uint8_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"u";
                        break;
                    }
                    case ArrayType::INT16: {
                        int16_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"i";
                        break;
                    }
                    case ArrayType::UINT16: {
                        uint16_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"u";
                        break;
                    }
                    case ArrayType::INT32: {
                        int32_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"i";
                        break;
                    }
                    case ArrayType::UINT32: {
                        uint32_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"u";
                        break;
                    }
                    case ArrayType::INT64: {
                        int64_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"i";
                        break;
                    }
                    case ArrayType::UINT64: {
                        uint64_t entry = col[iRow];
                        std::wstring wstr = std::to_wstring(entry);
                        rowStr += temp + std::u16string(wstr.begin(), wstr.end()) + u"u";
                        break;
                    }
                    case ArrayType::CELL: {
                        // cell array of char arrays
                        CharArrayRef entry = col[iRow][0];
                        String strEntry = entry.toUTF16();
                        rowStr += temp + u"\"" + static_cast<std::u16string>(strEntry) + u"\"";
                        break;

                    }
                    case ArrayType::MATLAB_STRING: {
                        String entry = col[iRow];
                        rowStr += temp + u"\"" + static_cast<std::u16string>(entry) + u"\"";
                        break;
                    }
                    case ArrayType::DOUBLE:
                    case ArrayType::SINGLE: {
                        double entry = col[iRow];
                        if (!std::isfinite(entry)) {
                            // entry is nan or inf
                            nanNum += 1;
                            is_nan = true;
                            // if all the fields in the row is of the value NaN, then not adding the entire row
                            if (fieldNames.end() - fieldNames.begin() == nanNum) addRow = false;
                        }
                        else {
                            std::wstring wstr = std::to_wstring(entry);
                            rowStr += temp + std::u16string(wstr.begin(), wstr.end());
                        }
                        break;
                    }
                    default: 
                    {
                        //handle the missing data type
                        nanNum += 1;
                        is_nan = true;
                        if (fieldNames.end() - fieldNames.begin() == nanNum) addRow = false;
                        break;
                    }
                }
                if (it == lastMinusOne && is_nan) {
                    // last field is a nan
                    rowStr.pop_back();
                    rowStr += u" ";
                }
                else if (it == lastMinusOne) rowStr += u" "; //last field but not nan
                else if (!is_nan) rowStr += u","; // not last field and not nan
            }

            if (addRow) {
                std::wstring wstr = std::to_wstring(timeDataArray[iRow]);
                res += rowStr + std::u16string(wstr.begin(), wstr.end()) + u"\n";
            }

        }

        ArrayFactory factory;
        TypedArray<MATLABString> outString = factory.createScalar(std::move(res));
        outputs[0] = outString;

    }
};

